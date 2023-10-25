# step_1.qc.R
# step 1: perform QC on cells
# 

library(Seurat)
library(magrittr)
library(future)

workdir <- "."
srcdir <- file.path(workdir, "src")
sourcedir <- file.path(workdir, "data")
umidir <- file.path(sourcedir, "UMI_counts")
figdir <- file.path(workdir, "result", "figure")
infodir <- file.path(workdir, "result", "info")

# project name
project <- "xiaohua"

# pattern for defining mitochondrial and ribosomal genes
mito.pattern <- "^mt-"
ribo.pattern <- "^Rpl|^Rps"

# pattern for Multiplexing Capture label
multi.pattern <- '^CMO'

# set a random seed
set.seed(98)

# load functions
source(file.path(srcdir, "my_functions.R"))

# set parallelization in Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=8*1024^3)

# sample information
sample.info <- data.frame(SeqName=c("WT","KC","KPC","KSC","KPSC","WT-P","KC-P","KPC-P","KSC-P","KPSC-P"), 
                          Name=c("WT","KC","KPC","KSC","KPSC","WT-P","KC-P","KPC-P","KSC-P","KPSC-P"), 
                          Mutation=rep(c('WT','KC','KPC','KSC','KPSC'),2), 
                          Treatment=c(rep('Untreated',5),rep('Treated',5)))
rownames(sample.info) <- c("D01","D02","D03","D04","D05","D06","D07","D08","D09","D10")

# load UMI counts matrix for each sample
raw.counts.list <- list()
for (k in 1:nrow(sample.info)){
  pid <- rownames(sample.info)[k]
  sid <- sample.info$SeqName[k]
  raw.counts.list[[k]] <- my.Read10X(file.path(sourcedir, sid), pid)
}
names(raw.counts.list) <- rownames(sample.info)

# merge UMI counts matrix
raw.counts.all <- my.MergeMatrix.v2(raw.counts.list)
# remove multiplexing capture labels
raw.counts.all <- raw.counts.all[grep(multi.pattern, rownames(raw.counts.all), invert=T, value=T),]

# Initialize the Seurat object with the raw (non-normalized data).
panc <- CreateSeuratObject(counts=raw.counts.all, project=project, assay="RNA", min.cells=0, min.features=0, 
                           names.field=1, names.delim="_", meta.data=NULL)

# calculates the mitochondrial/ribosomal genes per cell
panc[["percent.mito"]] <- PercentageFeatureSet(panc, pattern=mito.pattern)
panc[["percent.ribo"]] <- PercentageFeatureSet(panc, pattern=ribo.pattern)

# add sample information
tmeta <- data.frame(row.names=rownames(panc@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc@meta.data[,"orig.ident"])])
}
panc %<>% AddMetaData(metadata=tmeta)

# perform cell filtering
# nGene > 1000, nGene <= 6000, nUMI > 2000, nUMI <= 40000, percent.mito < 5%
panc %<>% subset(subset=nFeature_RNA > 1000 & nFeature_RNA <= 6000 & nCount_RNA > 2000 & 
                   nCount_RNA <= 40000 & percent.mito < 5)

# save Seurat object
saveRDS(panc, file=file.path(infodir, "panc.rds"))

# free space
rm(raw.counts.all)
rm(tmeta)
rm(tdic)
gc()
