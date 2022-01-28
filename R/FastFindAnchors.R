#' @import Matrix
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @export
FastFindAnchors = function(
  tmp.dir = NULL,
  nCores = NULL,
  k.filter = 100,
  k.anchor = 5
) {

  nCores = nCores %||% parallel::detectCores()

  object.ncells = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/others/object_ncells.rds"))
  offsets = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/others/offsets.rds"))
  features = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/others/features.rds"))
  nSample = length(object.ncells)

  combinations = expand.grid(1:nSample, 1:nSample)
  combinations = combinations[combinations$Var1 < combinations$Var2, , drop = FALSE]

  message("Finding all pairwise anchors")

  if (nSample > 100) {
    rna.list = pbmcapply::pbmclapply(
      1:nSample, function(i) {
        rna = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/", i,".rds"))
        rna = rna@assays$RNA@data[features,]
        return(rna)
      }
    )

    ncell = c()
    for (i in 1:length(rna.list)) {
      ncell = append(ncell, ncol(rna.list[[i]]))
    }

    rna.average = pbmcapply::pbmclapply(
      rna.list, function(rna) {
        rna = rowMeans(rna)
        return(rna)
      }, mc.cores = 20
    )
    rna.average = do.call(cbind, rna.average)
    rna.cor = cor(rna.average)

    sample.tree = Seurat:::BuildSampleTree(rna.cor)

    paired = c(1, 2)
    for (j in 1:nrow(sample.tree)) {
      idx = Seurat:::ParseMergePair(sample.tree, j)
      if (sum(ncell[idx$object2]) > sum(ncell[idx$object1])) {
        tmp = idx$object2
        idx$object2 = idx$object1
        idx$object1 = tmp
      }
      p = lapply(idx$object2, function(n) {
        a = idx$object1[which(rank(-rna.cor[n, idx$object1]) <= 5)]
        return(cbind(n, a))
      })
      paired = rbind(paired, do.call(rbind, p))
    }
    paired = paired[-1,]
    paired = data.frame(paired)
    colnames(paired) = c("Var1", "Var2")
    rm(rna.list)
    gc(reset = T)
    combinations = paired[!duplicated(paired),]
  }

  all.anchors = pbmcapply::pbmclapply(
    X = 1:nrow(x = combinations),
    FUN = function(row) {
      FindAnchorsPair(
        combinations = combinations,
        row = row,
        offsets = offsets,
        k.filter = k.filter,
        k.anchor = k.anchor,
        tmp.dir = tmp.dir,
        features = features
      )
    },
    mc.cores = nCores
  )

  all.anchors = do.call(what = 'rbind', args = all.anchors)
  all.anchors = rbind(all.anchors, all.anchors[, c(2, 1, 3)])
  # AddDatasetID = getFromNamespace("AddDatasetID", "Seurat")
  all.anchors = Seurat:::AddDatasetID(anchor.df = all.anchors, offsets = offsets, obj.lengths = object.ncells)
  all.anchors = data.table::as.data.table(all.anchors)
  data.table::setindex(all.anchors, dataset2)
  data.table::setindex(all.anchors, dataset1)

  message("Merging data")
  if (nSample <= 100) {
    anchor.group <- all.anchors %>% group_by(dataset1, dataset2) %>% summarise(n = n())
    similarity.matrix = matrix(data = 0, ncol = nSample, nrow = nSample)
    similarity.matrix[lower.tri(similarity.matrix, diag=FALSE) | upper.tri(similarity.matrix, diag=FALSE)] = anchor.group$n


    ncell.a = matrix(data = NA, ncol = nSample, nrow = nSample)
    for (i in 1:nSample) {
      ncell.a[i,] = object.ncells[i]
    }

    ncell.b = matrix(data = NA, ncol = nSample, nrow = nSample)
    for (i in 1:nSample) {
      ncell.b[,i] = object.ncells[i]
    }
    ncel.min = pmin(ncell.a, ncell.b)
    similarity.matrix = similarity.matrix/ncel.min
    similarity.matrix[upper.tri(x = similarity.matrix, diag = TRUE)] = NA
    # BuildSampleTree = getFromNamespace("BuildSampleTree", "Seurat")
    sample.tree = Seurat:::BuildSampleTree(similarity.matrix)

    rna.bind = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/1.rds"))
    rna.bind = rna.bind@assays$RNA@data[features,]
    for (i in 2:nSample) {
      rna = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/", i,".rds"))
      rna.bind = cbind(rna.bind, rna@assays$RNA@data[features,])
    }
  } else {
    n = ceiling(nSample/100)
    rna.list = pbmcapply::pbmclapply(
      1:n, function(i) {
        start = (i-1)*100 + 1
        end = min(nSample, i*100)
        rna.bind = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/", start,".rds"))
        rna.bind = rna.bind@assays$RNA@data[features,]
        for (j in (start+1):end) {
          rna = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/", j,".rds"))
          rna.bind = cbind(rna.bind, rna@assays$RNA@data[features,])
        }
        return(rna.bind)
      }, mc.cores = n
    )

    rna.bind = rna.list[[1]]
    for (i in 2:n) {
      rna.bind = cbind(rna.bind, rna.list[[i]])
    }
  }

  rna.bind = ScaleData(rna.bind, verbose = F)
  pca = RunPCA(rna.bind, verbose = F)

  saveRDS(pca, paste0(tmp.dir, "/FastIntegrationTmp/others/raw_pca.rds"), compress = F)
  saveRDS(sample.tree, paste0(tmp.dir, "/FastIntegrationTmp/others/sample_tree.rds"), compress = F)
  saveRDS(all.anchors, paste0(tmp.dir, "/FastIntegrationTmp/anchors/anchors.rds"), compress = F)
  return(NULL)
}

FindAnchorsPair = function(
  tmp.dir = NULL,
  combinations,
  row,
  offsets,
  k.filter,
  k.anchor,
  features
) {

  i = combinations[row, 1]
  j = combinations[row, 2]

  object.list.1 = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/", i, ".rds"))
  object.list.2 = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/", j, ".rds"))

  object.list.1@assays$RNA@counts = object.list.1@assays$RNA@data
  object.list.2@assays$RNA@counts = object.list.2@assays$RNA@data

  features = c(Seurat::VariableFeatures(object.list.1),Seurat::VariableFeatures(object.list.2))
  features = names(sort(-table(features)))[1:2000]
  object.1 = subset(x = object.list.1, features = features)
  object.2 = subset(x = object.list.2, features = features)

  ncell1 = ncol(object.1)
  ncell2 = ncol(object.2)

  object.1 = Seurat::ScaleData(object.1, features = features)
  object.2 = Seurat::ScaleData(object.2, features = features)

  object.1[["ToIntegrate"]] = object.1[["RNA"]]
  DefaultAssay(object = object.1) = "ToIntegrate"

  object.2[["ToIntegrate"]] = object.2[["RNA"]]
  DefaultAssay(object = object.2) = "ToIntegrate"


  object.pair = Seurat::RunCCA(
    object1 = object.1,
    object2 = object.2,
    assay1 = "ToIntegrate",
    assay2 = "ToIntegrate",
    features = features,
    renormalize = FALSE,
    rescale = FALSE,
    num.cc = 30
  )

  object.pair = Seurat::L2Dim(object = object.pair, reduction = "cca")

  anchors = Seurat:::FindAnchors(
    object.pair = object.pair,
    assay = c("ToIntegrate", "ToIntegrate"),
    slot = "data",
    cells1 = colnames(x = object.1),
    cells2 = colnames(x = object.2),
    internal.neighbors = list(NULL, NULL),
    reduction = "cca.l2",
    nn.reduction = "cca",
    k.filter = k.filter,
    k.anchor = k.anchor,
    nn.method = "rann",
    dims = 1:30
  )
  anchors[, 1] = anchors[, 1] + offsets[i]
  anchors[, 2] = anchors[, 2] + offsets[j]
  return(anchors)
}

#' @export
BuildIntegrationFile = function(
  rna.list,
  tmp.dir = NULL,
  nCores = NULL,
  nfeatures = 2000
) {
  message("Checking the rna list")
  rna.list = SeuratObject:::CheckDuplicateCellNames(rna.list)
  nCores = nCores %||% parallel::detectCores()

  if (is.null(tmp.dir)) {
    message("Not specify the temporary directory, will create one in current directory")
    tmp.dir = "./"
  }

  if(file.exists(tmp.dir) == F) {
    stop("The temporary directory not exist")
  } else {
    dir.create(paste0(tmp.dir, "/FastIntegrationTmp"))
    dir.create(paste0(tmp.dir, "/FastIntegrationTmp/raw"))
    dir.create(paste0(tmp.dir, "/FastIntegrationTmp/anchors"))
    dir.create(paste0(tmp.dir, "/FastIntegrationTmp/others"))
    dir.create(paste0(tmp.dir, "/FastIntegrationTmp/inte"))
  }

  message("Splitting and saving rna list in tmp dir")
  pbmcapply::pbmclapply(
    1:length(rna.list), function(i) {
      rna = rna.list[[i]]
      rna@assays$RNA@counts = matrix(0)
      saveRDS(rna, paste0(tmp.dir, "/FastIntegrationTmp/raw/", i, ".rds"), compress = F)
      return(NULL)
    }, mc.cores = nCores
  )

  features = Seurat::SelectIntegrationFeatures(rna.list, nfeatures = nfeatures, verbose = F)
  object.ncells = sapply(X = rna.list, FUN = function(x) dim(x = x)[2])
  offsets = as.vector(x = cumsum(x = c(0, object.ncells)))[1:length(x = rna.list)]

  saveRDS(features, paste0(tmp.dir, "/FastIntegrationTmp/others/features.rds"))
  saveRDS(object.ncells, paste0(tmp.dir, "/FastIntegrationTmp/others/object_ncells.rds"))
  saveRDS(offsets, paste0(tmp.dir, "/FastIntegrationTmp/others/offsets.rds"))
  return(NULL)
}

