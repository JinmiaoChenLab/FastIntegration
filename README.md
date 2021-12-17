# FastIntegration v1.0.0
FastIntegration is a fast and high-capacity version of Seurat Integration. FastIntegrate can integrate thousands of scRNA-seq datasets and outputs batch-corrected values for downstream analysis.


## Requirement
FastIntegration requires the following packages:

* [R](https://www.r-project.org/) (>= 4.0.0)
* [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) (>= 4.0.0)
* [SeuratObject](https://cran.r-project.org/web/packages/SeuratObject/index.html) (>= 4.0.0)
* [data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
* [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
* [tictoc](https://cran.r-project.org/web/packages/tictoc/index.html)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
* [pbmcapply](https://cran.r-project.org/web/packages/pbmcapply/index.html)

## Installation

```R
devtools::install_github("git@github.com:JinmiaoChenLab/FastIntegrate.git")
```

## Usage

### Onestop function
```R
library(FastIntegration)
# rna.list is the list of seurat object
data = OneStopIntegration(
  rna.list = rna.list, 
  tmp.dir = "../test/", 
  max.cores = 30,
  feature.to.integration = NULL,
  feature.to.return = NULL
)

```

### Step by step integration
```R
library(Seurat)
library(pbmcapply)
library(FastIntegration)

# rna.list is the list of seurat object
BuildIntegrationFile(rna.list = rna.list, tmp.dir = "./", nCores = 50)
FastFindAnchors(tmp.dir = "./", nCores = 50)

genes = readRDS("FastIntegrationTmp/raw/1.rds")
genes = rownames(genes)
idx = split(1:length(genes), cut(1:length(genes), 20, labels = FALSE))
pca = readRDS("raw_pca.rds")

pbmclapply(
  1:20, function(i) {
    rna.integrated = FastIntegration(tmp.dir = "./", npcs = 1:30, slot = "data",
                                     features.to.integrate = genes[idx[[i]]])
    saveRDS(rna.integrated, paste0("FastIntegrationTmp/inte/inte_", i, ".rds"), compress = F)
  }, mc.cores = 20
)

```


## Usage Scenario
We have apply FastIntegration to [DISCO](http://www.immunesinglecell.org/) database for integrating thousands of samples.

## License
All other code in this repository is licensed under a [GPL-3](https://www.r-project.org/Licenses/GPL-3) license.

