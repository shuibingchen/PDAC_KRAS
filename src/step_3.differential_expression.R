# step_3.differential_expression.R
# step 3: conduct differential expression analysis; prepare input files for pySCENIC
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

# set a random seed
set.seed(98)

# load functions
source(file.path(srcdir, "my_functions.R"))

# set parallelization in Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=8*1024^3)

# load Seurat object
panc <- readRDS(file.path(infodir, "panc.rds"))

# set cell groups
Idents(panc) <- 'Name'

# get number of cells per group
nCells <- FetchData(panc, vars=c('ident')) %>% dplyr::count(ident)

# conduct DE analysis between KC,KPC,KSC,KPSC sample and WT sample
tde <- my.DE.pair(panc, 'KC', 'WT', 20, nCells, 'KC.vs.WT')
tde <- my.DE.pair(panc, 'KPC', 'WT', 20, nCells, 'KPC.vs.WT')
tde <- my.DE.pair(panc, 'KSC', 'WT', 20, nCells, 'KSC.vs.WT')
tde <- my.DE.pair(panc, 'KPSC', 'WT', 20, nCells, 'KPSC.vs.WT')

# conduct DE analysis between treated and untreated sample pairs
tde <- my.DE.pair(panc, 'WT-P', 'WT', 20, nCells, 'WT-P.vs.WT-C')
tde <- my.DE.pair(panc, 'KC-P', 'KC', 20, nCells, 'KC-P.vs.KC-C')
tde <- my.DE.pair(panc, 'KPC-P', 'KPC', 20, nCells, 'KPC-P.vs.KPC-C')
tde <- my.DE.pair(panc, 'KSC-P', 'KSC', 20, nCells, 'KSC-P.vs.KSC-C')
tde <- my.DE.pair(panc, 'KPSC-P', 'KPSC', 20, nCells, 'KPSC-P.vs.KPSC-C')

# prepare input file for running pySCENIC
# apply a soft filtering on genes
my.genes <- myFilterGenes(tobj=panc, tassay="RNA", tclusts=NULL, tCountsPerGene=3*150, tminSamples=150)

# write filtered UMI counts matrix to file
write.table(as.data.frame(as.matrix(panc[["RNA"]]@counts[my.genes,])), 
            file=file.path(infodir,"UMI.counts.cleaned.csv"), quote=F, sep=",", row.names=T, col.names=NA)

