# step_2.cell_clustering.R
# step 2: perform data normalization and cell clustering
# 

library(Seurat)
library(magrittr)
library(future)
library(dplyr)
library(tibble)

workdir <- "."
srcdir <- file.path(workdir, "src")
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "result", "figure")
infodir <- file.path(workdir, "result", "info")

# pattern for defining mitochondrial and ribosomal genes
mito.pattern <- "^mt-"
ribo.pattern <- "^Rpl|^Rps"

# set a random seed
set.seed(98)

# load functions
source(file.path(srcdir, "my_functions.R"))

# set parallelization in Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=8*1024^3)

# load Seurat object
panc <- readRDS(file.path(infodir, "panc.rds"))

# normalization with sctransform
panc %<>% SCTransform(variable.features.n=3500, return.only.var.genes=TRUE)

# remove mitochondrial and ribosomal genes from variable gene list
# and select the remaining top 3000 genes for clustering analysis
variable.genes <- setdiff(VariableFeatures(panc), 
                          c(grep(mito.pattern, rownames(panc), value=T), grep(ribo.pattern, rownames(panc), value=T)))
VariableFeatures(panc) <- head(variable.genes, 3000)

# perform linear dimensional reduction
panc %<>% RunPCA(features=VariableFeatures(object=panc))

# select top 30 PCs
pcs <- 30

# perform cell clustering
panc %<>% FindNeighbors(reduction="pca", dims=1:pcs)
panc %<>% FindClusters(resolution=seq(0.1,1,by=0.1), verbose=T)

# perform dimensional reduction
panc %<>% RunUMAP(dims=1:pcs, reduction="pca", n.components=3, seed.use=42, n.neighbors=30, n.epochs=1000)

# set cell identity
panc %<>% SetIdent(value="SCT_snn_res.0.1")

# reorder clusters
Idents(panc) <- factor(Idents(panc), levels=0:(length(unique(Idents(panc)))-1))

# save Seurat object
saveRDS(panc, file=file.path(infodir, 'panc.rds'))
