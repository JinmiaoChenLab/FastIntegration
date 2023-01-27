GetDiscoSample = function(
    verbose = T
) {
  if(!require("jsonlite")) {
    stop("Please install jsonlite")
  }

  tryCatch(
    {
      if (verbose) {
        message("Starting to download metadata from DISCO database")
      }
      meta = fromJSON("https://www.immunesinglecell.org/api/vishuo/sample/all")
      return(meta)
    },
    error=function(cond) {
      stop("Network error. Please try again")
    }
  )
}

#' @import jsonlite
#' @export
FindSampleByMetadata = function(
    tissue = c(),
    disease = c(),
    platform = c(),
    project.id = c(),
    sample.id = c(),
    sample.type = c()
) {
  meta.all =  GetDiscoSample()
  if (length(tissue) > 0) {
    meta.all = meta.all[which(meta.all$tissue %in% tissue),]
  }
  if (length(disease) > 0) {
    meta.all = meta.all[which(meta.all$disease %in% disease),]
  }
  if (length(platform) > 0) {
    meta.all = meta.all[which(meta.all$platform %in% platform),]
  }
  if (length(project.id) > 0) {
    meta.all = meta.all[which(meta.all$projectId %in% project.id),]
  }
  if (length(sample.type) > 0) {
    meta.all = meta.all[which(meta.all$sampleType %in% sample.type),]
  }

  meta.all = meta.all[which(meta.all$processStatus == "QC pass"),]
  if (nrow(meta.all) == 0) {
    stop("Sorry, no sample is found. Please try to use other filters.")
  } else {
    return(meta.all)
  }
}

#' @export
DownloadDiscoData = function(
    metadata,
    dir = "./disco_data"
) {
  options(timeout = 20*60)
  meta = metadata
  dir.create(dir)
  tryCatch(
    {
      for (i in 1:nrow(meta)) {
        outpur.file = paste0(dir, "/", meta$sampleId[i], ".rds")
        if (file.exists(outpur.file) & tools::md5sum(outpur.file) == meta$md5[i]) {
          message(paste0(i, "st sample has been downloaded before. Ignore..."))
        } else {
          message(paste0("Downloading the ", i, "st sample"))
          # rna = readRDS(url(paste0(
          #   "http://dc.vishuo.com:8887/api/vishuo/download/getExp?project=",meta$projectId[i],"&sample=",meta$sampleId[i]
          # )))
          download.file(paste0("http://dc.vishuo.com:8887/api/vishuo/download/getExp?project=",meta$projectId[i],"&sample=",meta$sampleId[i]),
                        paste0(dir, "/", meta$sampleId[i], ".rds"))
          # saveRDS(rna, paste0(dir, "/", meta$sampleId[i], ".rds"), compress = F)
        }
      }
    },
    error=function(cond) {
      stop("Network error. Please try again")
    }
  )
  message("Job finished")
}


AddCountsSlot = function(rna) {
  rna@assays$RNA@counts = sweep(expm1(rna@assays$RNA@data), 2, rna$nCount_RNA, `*`) /10000
  return(rna)
}



