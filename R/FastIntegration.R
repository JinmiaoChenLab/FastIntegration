#' @import tictoc
#' @importFrom dplyr `%>%`
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @import RcppEigen
#' @import Rcpp
#' @rawNamespace useDynLib(FastIntegration)
#' @export
FastIntegration = function(
  tmp.dir = NULL,
  features.to.integrate = NULL,
  npcs = 1:30,
  nn.k = 100,
  slot = c("data", "counts")
) {
  data.table::setDTthreads(threads = 1L)

  input.pca = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/others/raw_pca.rds"))
  input.pca = input.pca@cell.embeddings[,npcs]
  message("Reading anchor file")
  anchors = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/anchors/anchors.rds"))
  offsets = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/others/offsets.rds"))
  obj.lengths = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/others/object_ncells.rds"))
  sample.tree = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/others/sample_tree.rds"))

  message("Reading rna list file")
  object.list = list()
  for (i in 1:length(obj.lengths)) {
    a = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/", i, ".rds"))
    if (slot == "data") {
      object.list[[i]] = a@assays$RNA@data[features.to.integrate,]
    }else {
      object.list[[i]] = exp(a@assays$RNA@data[features.to.integrate,])-1
      object.list[[i]] = drop0(object.list[[i]])
    }
  }

  id.table = cbind(paste0("1_", 1:ncol(object.list[[1]])), colnames(object.list[[1]]))
  for (i in 2:length(object.list)){
    id.table = rbind(id.table, cbind(paste0(i, "_", 1:ncol(object.list[[i]])), colnames(object.list[[i]])))
  }
  id.table = data.frame(id.table)
  rownames(id.table) = id.table$X1

  reference.datasets = 1:length(object.list)

  reference.integrated = FastPairwiseIntegrateReference(
    anchors = anchors,
    offsets = offsets,
    reference.objects = reference.datasets,
    object.list = object.list,
    features.to.integrate = features.to.integrate,
    input.pca = input.pca,
    id.table = id.table,
    nn.k = nn.k,
    sample.tree = sample.tree,
    objects.ncell = obj.lengths
  )
}

FastPairwiseIntegrateReference = function(
  anchors,
  offsets,
  reference.objects,
  object.list,
  features.to.integrate = NULL,
  input.pca ,
  id.table,
  nn.k,
  sample.tree,
  objects.ncell
) {

  cellnames.list = list()
  for (ii in 1:length(object.list)) {
    cellnames.list[[ii]] = colnames(x = object.list[[ii]])
  }
  names(object.list) = as.character(-(1:length(object.list)))

  for (ii in 1:nrow(sample.tree)) {
    merge.pair = as.character(sample.tree[ii, ])
    length1 = ncol(object.list[[merge.pair[1]]])
    length2 = ncol(object.list[[merge.pair[2]]])

    if (length2 > length1) {
      merge.pair = rev(x = merge.pair)
      sample.tree[ii, ] = as.numeric(merge.pair)
    }
    ParseMergePair = getFromNamespace("ParseMergePair", "Seurat")
    datasets = ParseMergePair(sample.tree, ii)

    object.1 = object.list[[merge.pair[1]]]
    object.2 = object.list[[merge.pair[2]]]

    tictoc::tic("Filtering anchors")
    filtered.anchors = data.frame(anchors[dataset1 %in% datasets$object1 & dataset2 %in% datasets$object2,])
    tictoc::toc()

    message(
      "Merging dataset ",
      paste(datasets$object2, collapse = " "),
      " into ",
      paste(datasets$object1, collapse = " ")
    )
    message(paste0("length1:", length1, "," , "length2:", length2))

    integrated.matrix = FastRunIntegration(
      filtered.anchors = filtered.anchors,
      reference = object.1,
      query = object.2,
      cellnames.list = cellnames.list,
      features.to.integrate = features.to.integrate,
      input.pca = input.pca,
      id.table = id.table,
      nn.k = nn.k
    )

    object.list[[as.character(ii)]] = integrated.matrix
    object.list[[merge.pair[[1]]]] = NULL
    object.list[[merge.pair[[2]]]] = NULL
  }

  return(object.list[[as.character(ii)]])
}

FastRunIntegration = function(
  filtered.anchors,
  reference,
  query,
  cellnames.list,
  features.to.integrate,
  input.pca,
  id.table,
  nn.k
) {

  if (dim(filtered.anchors)[1] > 100000) {
    filtered.anchors = filtered.anchors %>%
      arrange(desc(score)) %>%
      group_by(cell2) %>% slice(1:10)
    filtered.anchors = data.frame(filtered.anchors)
  }

  tictoc::tic("Getting cell offsets")
  cells1 = colnames(x = reference)
  cells2 = colnames(x = query)
  cell1.offset = match(id.table[paste0(filtered.anchors$dataset1, "_", filtered.anchors$cell1),"X2"], cells1)
  cell2.offset = match(id.table[paste0(filtered.anchors$dataset2, "_", filtered.anchors$cell2),"X2"], cells2)
  filtered.anchors[, 1] = cell1.offset
  filtered.anchors[, 2] = cell2.offset
  tictoc::toc()

  tictoc::tic("Weighting")
  weight.matrix = FastFindWeights(
    cells1,
    cells2,
    anchors = filtered.anchors,
    reduction = input.pca,
    nn.k = nn.k
  )
  tictoc::toc()

  anchors1 = cells1[filtered.anchors[, "cell1"]]
  anchors2 = cells2[filtered.anchors[, "cell2"]]

  tictoc::tic("integrating data")
  # IntegrateDataC = getFromNamespace("IntegrateDataC", "Seurat")
  integration.matrix = query[, anchors2] - reference[, anchors1]
  if (dim(query)[2] > 100000) {
    ss = split(1:dim(query)[2], cut(seq_along(1:dim(query)[2]), ceiling(dim(query)[2]/20000), labels = FALSE))
    integrated = pbmcapply::pbmclapply(
      1:length(ss),
      function(i) {
        integrated = query[,ss[[i]]] - integration.matrix %*% weight.matrix[,ss[[i]]]
        # integrated = IntegrateDataC(integration.matrix, weight.matrix[,ss[[i]]], query[,ss[[i]]])
      },
      mc.cores = min(ceiling(dim(query)[2]/50000), 40)
    )
    integrated = do.call("cbind", integrated)
  }else {
    integrated = query - integration.matrix %*% weight.matrix
    # integrated = IntegrateDataC(integration.matrix, weight.matrix, query)
  }
  # integrated[which(integrated < 0.2)] = 0
  # integrated = Matrix::drop0(integrated)
  dimnames(integrated) = dimnames(query)
  new.expression = cbind(reference, integrated)

  tictoc::toc()

  return(new.expression)
}


FastFindWeights = function(
  nn.cells1,
  nn.cells2,
  anchors,
  reduction = NULL,
  nn.k = nn.k
) {

  anchors.cells2 = unique(x = nn.cells2[anchors[, "cell2"]])
  data.use = reduction[nn.cells2, ]

  NNHelper = getFromNamespace("NNHelper", "Seurat")
  nn.k = min(nn.k, length(anchors.cells2))
  knn_2_2 = NNHelper(
    data = data.use[anchors.cells2, ],
    query = data.use,
    k = nn.k,
    method = "rann",
    n.trees = 50,
    eps = 0
  )

  distances = Seurat::Distances(object = knn_2_2)
  distances = 1 - (distances / distances[, ncol(x = distances)])
  cell.index = Seurat::Indices(object = knn_2_2)
  FindWeightsC = getFromNamespace("FindWeightsC", "Seurat")
  weights = FindWeightsC(
    cells2 = 0:(length(x = nn.cells2) - 1),
    distances = as.matrix(x = distances),
    anchor_cells2 = anchors.cells2,
    integration_matrix_rownames = nn.cells2[anchors$cell2],
    cell_index = cell.index,
    anchor_score = anchors[, "score"],
    min_dist = 0,
    sd = 1,
    display_progress = TRUE
  )

  return(weights)
}

