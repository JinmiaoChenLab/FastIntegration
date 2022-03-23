#' @export
OneStopIntegration = function(
  rna.list = NULL,
  tmp.dir = "./",
  feature.to.integration = NULL,
  feature.to.return = NULL,
  max.cores = 30,
  nfeatures = 2000
) {
  nCores = max.cores
  message("Building integration files")
  BuildIntegrationFile(rna.list = rna.list, tmp.dir = tmp.dir, nCores = nCores, nfeatures = nfeatures)

  message("Finding anchors")
  FastFindAnchors(tmp.dir = tmp.dir, nCores = nCores)

  message("Integrating data")
  if (is.null(feature.to.integration) == T) {
    feature.to.integration = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/raw/1.rds"))
    feature.to.integration = rownames(feature.to.integration)
    idx = split(1:length(feature.to.integration), cut(1:length(feature.to.integration), 20, labels = FALSE))
  }

  pbmcapply::pbmclapply(
    1:20, function(i) {
      rna.integrated = FastIntegration(tmp.dir = tmp.dir, npcs = 1:30, slot = "data",
                                       features.to.integrate = feature.to.integration[idx[[i]]])
      saveRDS(rna.integrated, paste0(tmp.dir, "/FastIntegrationTmp/inte/inte_", i, ".rds"), compress = F)
    }, mc.cores = 20
  )

  features = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/others/features.rds"))
  feature.to.return = feature.to.return %||% features

  rna.data = pbmcapply::pbmclapply(
    1:20, function(i) {
      rna = readRDS(paste0(tmp.dir, "/FastIntegrationTmp/inte/inte_", i, ".rds"))
      if (length(intersect(rownames(rna), feature.to.return)) > 0) {
        rna = rna[intersect(feature.to.return, rownames(rna)), , drop = F]
        return(rna)
      } else {
        return(NULL)
      }
    }, mc.cores = 20
  )
  na.id = unlist(lapply(rna.data, is.null))
  rna.data = do.call(rbind, rna.data[which(na.id == F)])
  rna.data = CreateSeuratObject(rna.data)
  return(rna.data)
}
