#' @import Seurat
#' @import SeuratObject
#' @import data.table
#' @import dplyr
#' @import Matrix
#' @importFrom utils getFromNamespace
#' @export
FastFindAnchors = function(
  tmp.dir = NULL,
  nCores = NULL,
  k.filter = 200,
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
  AddDatasetID = getFromNamespace("AddDatasetID", "Seurat")
  all.anchors = AddDatasetID(anchor.df = all.anchors, offsets = offsets, obj.lengths = object.ncells)
  all.anchors = as.data.table(all.anchors)
  setindex(all.anchors, dataset2)
  setindex(all.anchors, dataset1)

  message("Splitting and saving rna list in tmp dir")
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
  BuildSampleTree = getFromNamespace("BuildSampleTree", "Seurat")
  sample.tree = BuildSampleTree(similarity.matrix)

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

  object.1 = subset(x = object.list.1, features = features)
  object.2 = subset(x = object.list.2, features = features)

  ncell1 = ncol(object.1)
  ncell2 = ncol(object.2)

  object.1 = ScaleData(object.1, features = features)
  object.2 = ScaleData(object.2, features = features)

  object.1[["ToIntegrate"]] = object.1[["RNA"]]
  DefaultAssay(object = object.1) = "ToIntegrate"

  object.2[["ToIntegrate"]] = object.2[["RNA"]]
  DefaultAssay(object = object.2) = "ToIntegrate"

  object.pair = RunCCA(
    object1 = object.1,
    object2 = object.2,
    assay1 = "ToIntegrate",
    assay2 = "ToIntegrate",
    features = features,
    renormalize = FALSE,
    rescale = FALSE
  )

  object.pair = L2Dim(object = object.pair, reduction = "cca")

  FindAnchors = getFromNamespace("FindAnchors", "Seurat")
  anchors = FindAnchors(
    object.pair = object.pair,
    assay = c("ToIntegrate", "ToIntegrate"),
    slot = "data",
    cells1 = colnames(x = object.1),
    cells2 = colnames(x = object.2),
    internal.neighbors = list(NULL, NULL),
    reduction = "cca.l2",
    nn.reduction = "cca",
    k.filter = k.filter,
    k.anchor = k.anchor
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
  CheckDuplicateCellNames = getFromNamespace("CheckDuplicateCellNames", "SeuratObject")
  rna.list = CheckDuplicateCellNames(rna.list)
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

  features = Seurat::SelectIntegrationFeatures(rna.list, nfeatures = nfeatures)
  object.ncells = sapply(X = rna.list, FUN = function(x) dim(x = x)[2])
  offsets = as.vector(x = cumsum(x = c(0, object.ncells)))[1:length(x = rna.list)]

  saveRDS(features, paste0(tmp.dir, "/FastIntegrationTmp/others/features.rds"))
  saveRDS(object.ncells, paste0(tmp.dir, "/FastIntegrationTmp/others/object_ncells.rds"))
  saveRDS(offsets, paste0(tmp.dir, "/FastIntegrationTmp/others/offsets.rds"))
  return(NULL)
}

