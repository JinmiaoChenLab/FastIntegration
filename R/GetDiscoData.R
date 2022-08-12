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
    expressed.gene = c(),
    unexpressed.gene = c(),
    dir = "./disco_data"
) {

  meta = metadata
  dir.create(dir)
  tryCatch(
    {
      for (i in 1:nrow(meta)) {
        message(paste0("Downloading the ", i, "st sample"))
        rna = readRDS(url(paste0(
          "http://dc.vishuo.com:8887/api/vishuo/download/getExp?project=",meta$projectId[i],"&sample=",meta$sampleId[i]
        )))

        if (length(expressed.gene) > 0) {
          for (j in 1:length(expressed.gene)) {
            if (length(which(rna@assays$RNA@data[expressed.gene[j],] > 0)) > 0) {
              rna = subset(rna, cells = names(which(rna@assays$RNA@data[expressed.gene[j],] > 0)))
            } else {
              rna = NULL
              break
            }
          }
        }

        if (length(unexpressed.gene) > 0) {
          for (j in 1:length(unexpressed.gene)) {
            if (length(which(rna@assays$RNA@data[unexpressed.gene[j],] == 0)) > 0) {
              rna = subset(rna, cells = names(which(rna@assays$RNA@data[unexpressed.gene[j],] == 0)))
            } else {
              rna = NULL
              break
            }
          }
        }

        if (is.null(rna) == F) {
          rna@assays$RNA@counts = expm1(rna@assays$RNA@data)
          saveRDS(rna, paste0(dir, "/", meta$sampleId[i], ".rds"), compress = F)
        } else {
          message(paste0("Sking the ", i, "st sample. No cells are found after filtering."))
        }

      }
    },
    error=function(cond) {
      stop("Network error. Please try again")
    }
  )
  message("Job finished")
}





