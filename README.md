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

We highly recommend you to build R with openblas which will accelerate integration 2-3x times.

Here is the common way to do it:

sudo yum install -y openblas openblas-threads openblas-openmp # for centos

sudo apt-get install libopenblas-dev # for debian

./configure --enable-R-shlib --enable-byte-compiled-packages \

  --enable-BLAS-shlib --enable-memory-profiling \
  
  --with-blas="-lopenblas"
  
## Installation

```R
devtools::install_github("git@github.com:JinmiaoChenLab/FastIntegrate.git")
```

## Usage

### Preprocess
```R
library(Seurat)
library(pbmcapply)
rna.list = readRDS("rna_list.rds") # read list of Seurat object, each element in list is a sample

# make all samples have same genes
overlapped.gene = Reduce(intersect, lapply(rna.list, rownames))
for (i in 1:length(rna.list)) {
  rna.list[[i]] = subset(rna.list[[i]], features = overlapped.gene)
  rna.list[[i]] = NormalizeData(rna.list[[i]])
  rna.list[[i]] = FindVariableFeatures(rna.list[[i]])
  rna.list[[i]] = RenameCells(rna.list[[i]], new.names = paste0(Cells(rna.list[[i]]), "--", i))
}

```



### Onestop function

For large sample size (> 200 samples), we recommend to use step by step integration.
```R
library(FastIntegration)
# rna.list is the list of seurat object
data = OneStopIntegration(
  rna.list = rna.list, 
  tmp.dir = "./test/", 
  max.cores = 30
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

# close current R session and open a new one to clean the memory (This is important for large data integration)
# In the new session, please just set work directory and do not load any data. Then run the following codes:

library(Seurat)
library(pbmcapply)
library(FastIntegration)

genes = readRDS("FastIntegrationTmp/raw/1.rds")
genes = rownames(genes)
idx = split(1:length(genes), cut(1:length(genes), 20, labels = FALSE))
pbmclapply(
  1:20, function(i) {
    rna.integrated = FastIntegration(tmp.dir = "./", npcs = 1:30, slot = "data",
                                     features.to.integrate = genes[idx[[i]]])
    saveRDS(rna.integrated, paste0("FastIntegrationTmp/inte/inte_", i, ".rds"), compress = F)
  }, mc.cores = 20
)

```


### After integration
```R
##### create Seurat obj with the variable features of integration (For very big dataset) ##### 
features = readRDS("FastIntegrationTmp/others/features.rds")
rna.data = pbmclapply(
  1:20, function(i) {
    rna = readRDS(paste0("./FastIntegrationTmp/inte/inte_", i, ".rds"))
    rna = rna[intersect(rownames(rna), features),]
    return(rna)
  }, mc.cores = 20
)
rna.data = do.call(rbind, rna.data)
rna.data = CreateSeuratObject(rna.data)
rna.data = ScaleData(rna.data, features = features)
rna.data = RunPCA(rna.data, features = features, npcs = 50)
rna.data = FindNeighbors(rna.data, dims = 1:50)
rna.data = FindClusters(rna.data, graph.name = "RNA_snn", algorithm = 2)
rna.data = RunUMAP(rna.data, dims = 1:50)


##### select varibale gene based on integrated data  (For dataset with less than 100 samples) #####
rna.data = pbmclapply(
  1:20, function(i) {
    rna = readRDS(paste0("FastIntegrationTmp/inte/inte_", i, ".rds"))
    return(rna)
  }, mc.cores = 20
)

rna.data = do.call(rbind, rna.data)
rna.data = CreateSeuratObject(rna.data)
rna.data = FindVariableFeatures(rna.data, nfeatures = 2000)
features = VariableFeatures(rna.data)
rna.data = ScaleData(rna.data, features = features)
rna.data = RunPCA(rna.data, features = features)
rna.data = FindNeighbors(rna.data, dims = 1:50)
rna.data = FindClusters(rna.data, resolution = 0.5, algorithm = 2)
rna.data = RunUMAP(rna.data, dims = 1:50)

```


## Usage Scenario
We have apply FastIntegration to [DISCO](http://www.immunesinglecell.org/) database for integrating thousands of samples.

## License
All other code in this repository is licensed under a [GPL-3](https://www.r-project.org/Licenses/GPL-3) license.

