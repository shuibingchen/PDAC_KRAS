# my_functions.R
# customized functions for processing data and plotting
# 

library(Seurat, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(Matrix, quietly=T, warn.conflicts=F)
library(scater, quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pheatmap, quietly=T, warn.conflicts=F)
library(magrittr, quietly=T, warn.conflicts=F)
library(cowplot, quietly=T, warn.conflicts=F)
library(R.utils, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)
library(tibble, quietly=T, warn.conflicts=F)
library(ggrepel, quietly=T, warn.conflicts=F)
library(svglite, quietly=T, warn.conflicts=F)

# read in raw counts data from a given sample
my.Read10X <- function(tcountdir, tprefix=NULL, tmapfile=NULL){
  # directory exists?
  if (! dir.exists(tcountdir)){
    print(paste("input raw counts folder does NOT exist:", tcountdir, sep=" "))
    return(NULL)
  }
  # file exists?
  tmat.file <- paste(tcountdir, "matrix.mtx.gz", sep="/")
  tfnames.file <- paste(tcountdir, "features.tsv.gz", sep="/")
  tbnames.file <- paste(tcountdir, "barcodes.tsv.gz", sep="/")
  for (tf in c(tmat.file, tfnames.file, tbnames.file)){
    if (! file.exists(tf)){
      print(paste("input file does NOT exist:", tf))
      return(NULL)
    }
  }
  # extract counts matrix
  cat(paste("Loading UMI counts table from", tcountdir, "..."))
  tmat <- readMM(paste(tcountdir, "matrix.mtx.gz", sep="/"))
  tfnames <- read.delim(paste(tcountdir, "features.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  tbnames <- read.delim(paste(tcountdir, "barcodes.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  # update column names (cell ids)
  if (is.null(tprefix)){
    colnames(tmat) = tbnames$V1
  }
  else{
    colnames(tmat) = paste(tprefix, tbnames$V1, sep="_")
  }
  # replace rowname (Ensembl id) by gene symbol
  # in case gene symbol is not unique, append the _EnsemblID after it
  # missing gene symbol will be replaced by EnsemblID
  rownames(tmat) <- uniquifyFeatureNames(ID=tfnames$V1, names=tfnames$V2)
  if (! is.null(tmapfile)){
    tmap <- data.frame(ensembl_id=as.vector(tfnames$V1), name=as.vector(tfnames$V2), unique_name=rownames(tmat))
    write.table(tmap, file=tmapfile, quote=F, sep="\t", row.names=F, col.names=T)
  }
  cat(" done.","\n")
  return(tmat)
}

# merge raw read counts table collected from multiple samples, in a more efficient way
my.MergeMatrix.v2 <- function(tmats){
  cat("Merge raw UMI counts ")
  tfunc <- function(x,y){
    tres <- cbind(x,y[rownames(x),])
    cat(".")
    return(tres)
  }
  tmerged <- Reduce(f=tfunc, x=tmats)
  # fill na with 0
  tmerged[is.na(tmerged)] <- 0
  cat(" done.")
  return(tmerged)
}

# function to compare pre-defined two groups of cells
my.DE.pair <- function(tobj, cdtA, cdtB, ntop, nCells, ttitle, mincells=10, 
                       min.pct=0.1, logfc.threshold=0.1, test.use='wilcox'){
  # check number of cells in each condition
  nA <- as.numeric(nCells %>% filter(ident==cdtA) %>% select(n))
  nB <- as.numeric(nCells %>% filter(ident==cdtB) %>% select(n))
  print(paste0('Compare ',cdtA,' (',nA,' cells)',' to ',cdtB, ' (',nB,' cells).'))
  # enough cells?
  if (nA < mincells | nB < mincells){
    print('Too few cells to perform DE.')
    return(NULL)
  }
  # run DE
  de.results <- FindMarkers(tobj, ident.1=cdtA, ident.2=cdtB, min.pct=min.pct, logfc.threshold=logfc.threshold, test.use=test.use)
  # write to file
  out.file.prefix <- paste('DE',cdtA,'vs',cdtB,paste0('mincells_',mincells),
                           test.use,paste0('min_pct_',min.pct),paste0('logfc_',logfc.threshold), sep='.')
  de.file <- file.path(infodir, paste(out.file.prefix, 'txt', sep='.'))
  write.table(de.results, de.file, quote=F, sep='\t', col.names=NA)
  return(de.results)
}

# plot cells colored by an annotation (DimPlot)
myDimPlot <- function(tobj, treduct, tcate, torder=NULL, tsuffix, tcells=NULL, tcolor=NULL, 
                      tsort=FALSE, tlabel=FALSE, tsplit=FALSE, txlim=NULL, tylim=NULL,
                      tncol=4, tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20){
  tdataToPlot <- data.frame()
  tlabel.pos <- data.frame()
  tg <- ggplot()
  if (treduct == "tsne"){
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("tSNE_1","tSNE_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      missed.categroy <- setdiff(unique(tdataToPlot$Category), torder)
      if (length(missed.categroy) > 0){
        tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=missed.categroy))
      }
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    if (tsort){
      tdataToPlot <- tdataToPlot[with(tdataToPlot, order(Category)),]
    }
    tg <- ggplot(tdataToPlot, aes(x=tSNE_1, y=tSNE_2, color=Category))
    tg <- tg + ggtitle(paste("tSNE","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(tSNE_1, tSNE_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  } else if (treduct == "umap") {
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("UMAP_1","UMAP_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      missed.categroy <- setdiff(unique(tdataToPlot$Category), torder)
      if (length(missed.categroy) > 0){
        tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=missed.categroy))
      }
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    if (tsort){
      tdataToPlot <- tdataToPlot[with(tdataToPlot, order(Category)),]
    }
    tg <- ggplot(tdataToPlot, aes(x=UMAP_1, y=UMAP_2, color=Category))
    tg <- tg + ggtitle(paste("UMAP","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(UMAP_1, UMAP_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  }
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (! is.null(tcolor)){
    tg <- tg + scale_color_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (tsplit == TRUE){
    tg <- tg + facet_wrap(~Category, ncol=tncol)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Category), color="black")
  }
  return(tg)
}

# plot cells colored by an annotation (DimPlot)
# similar to myDimPlot, except that split by another label (e.g. color by ident but split by sample)
# this version is a more generalized function 'myDimPlot'
myDimPlot3 <- function(tobj, treduct, tgroup_by, tgroup_order=NULL, thighlight=NULL, tsuffix, 
                       tcells=NULL, tcolor=NULL, tlabel=FALSE, tsplit_by=NULL, tsplit_order=NULL, txlim=NULL, tylim=NULL,
                       tncol=1, tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20, tlbsize=2){
  # set coordinates variable name
  vars.reduct <- c("UMAP_1","UMAP_2")
  if (treduct == "tsne"){
    vars.reduct <- c("tSNE_1","tSNE_2")
  }
  # extract coordinates + group
  tdataToPlot <- FetchData(tobj, cells=tcells, vars=c(vars.reduct, tgroup_by))
  colnames(tdataToPlot) <- c("Dim_1","Dim_2","Group")
  # update group order if available
  if (!is.null(tgroup_order)){
    # add fake rows to make sure each group is considered
    tmiss.cates <- setdiff(tgroup_order, unique(tdataToPlot$Group))
    if (length(tmiss.cates) > 0){
      tdataToPlot <- rbind(tdataToPlot, data.frame(Dim_1=NA, Dim_2=NA, Group=tmiss.cates))
    }
    # reorder categories
    tdataToPlot$Group <- factor(tdataToPlot$Group, levels=tgroup_order)
  }
  # extract split
  if (! is.null(tsplit_by)){
    tsp <- FetchData(tobj, cells=tcells, vars=c(tsplit_by))
    tdataToPlot$Split <- tsp[rownames(tdataToPlot), c(tsplit_by)]
    # update split order if available
    if (! is.null(tsplit_order)){
      tdataToPlot$Split <- factor(tdataToPlot$Split, levels=tsplit_order)
    }
  }
  # reorder cells that needs to highlight (draw those cell points later)
  if (!is.null(thighlight)){
    tdataToPlot <- rbind(subset(tdataToPlot, ! Group %in% thighlight), subset(tdataToPlot, Group %in% thighlight))
  }
  # prepare group labeling
  tlabel.pos <- aggregate(cbind(Dim_1, Dim_2) ~ Group, data=tdataToPlot, FUN=median)
  colnames(tlabel.pos) <- c("Group","X","Y")
  # plot
  tg <- ggplot(tdataToPlot, aes(x=Dim_1, y=Dim_2, color=Group))
  tg <- tg + ggtitle(paste(toupper(treduct),"plots","by",tsuffix, sep=" "))
  # set range on coordinates
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (! is.null(tcolor)){
    tg <- tg + scale_color_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (! is.null(tsplit_by)){
    tg <- tg + facet_wrap(~Split, ncol=tncol)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Group), color="black", size=tlbsize)
  }
  return(tg)
}

# highlight expression of a set of genes (FeaturePlot)
MyFeaturePlot <- function(tobj, tgenes, tcells=NULL, tassay="RNA", treduction.name="umap", 
                          tsort=FALSE, txlim=NULL, tylim=NULL, tbreaks=NULL, tlimits=NULL, 
                          tlowcolor='gray80', thighcolor='red2', tncol=2, tlegend=NULL, tptsize=2, 
                          talpha=0.7, twidth=15, theight=12.5, tunits="in", tres=300){
  # genes valid?
  tgenes.valid <- intersect(tgenes, rownames(tobj[[tassay]]@data))
  if (length(tgenes.valid) == 0){
    cat("No valid genes found, do nothing!")
    return(NULL)
  }
  # assay valid?
  if (! tassay %in% names(tobj)){
    cat(paste("Not a valid assay:",tassay,sep=" "))
    return(NULL)
  }
  # extract gene expression
  texp <- as.matrix(tobj[[tassay]]@data[tgenes.valid, ,drop=F])
  # get coordinates
  tvars <- c("UMAP_1","UMAP_2")
  if (treduction.name == "tsne"){
    tvars <- c("tSNE_1","tSNE_2")
  }
  tdata <- FetchData(object=tobj, vars=tvars)
  colnames(tdata) <- c("X","Y")
  # plot
  tplots <- list()
  tk <- 1
  for (tgene in tgenes){
    # merge data for plotting
    tdata.merged <- merge(tdata, t(texp[tgene,,drop=F]), by=0, all=T)
    rownames(tdata.merged) <- tdata.merged$Row.names
    tdata.merged <- tdata.merged[,-1]
    colnames(tdata.merged) <- c("X","Y","Expression")
    # subset cells?
    if (! is.null(tcells)){
      tdata.merged <- tdata.merged[tcells,,drop=F]
    }
    # reorder cells by expression of the given gene
    if (tsort){
      tdata.merged <- tdata.merged[with(tdata.merged, order(Expression)),]
    }
    # plot
    if (max(tdata.merged$Expression) > 0){ # expressed in at least one cell
      # plot (rename x and y axis)
      tg <- ggplot(tdata.merged, aes(x=X, y=Y, color=Expression))
      tg <- tg + geom_point(shape=19, size=tptsize, alpha=talpha)
      if (! is.null(tbreaks)){
        tg <- tg + scale_color_gradient(low=tlowcolor, high=thighcolor, breaks=tbreaks, limits=tlimits)
      } else {
        tg <- tg + scale_color_gradient(low=tlowcolor, high=thighcolor)
      }
      if(! is.null(txlim)){
        tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
      }
      tg <- tg + ggtitle(tgene)
      if (treduction.name == "tsne"){
        tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
      } else {
        tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
      }
      tg <- tg + theme_bw()
      tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
      tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
      # add to list
      tplots[[tk]] <- tg
      tk <- tk + 1
    } else { # no expressions at all
      tg <- ggplot(tdata.merged, aes(x=X, y=Y))
      tg <- tg + geom_point(color="gray80", shape=19, size=tptsize, alpha=talpha)
      if(! is.null(txlim)){
        tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
      }
      tg <- tg + ggtitle(tgene)
      if (treduction.name == "tsne"){
        tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
      } else {
        tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
      }
      tg <- tg + theme_bw()
      tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
      tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
      # add to list
      tplots[[tk]] <- tg
      tk <- tk + 1
    }
  }
  # combine plots with Seurat::CombinePlots
  tcombined <- CombinePlots(tplots, ncol=tncol, legend=tlegend)
  return(tcombined)
}

# modify the DotPlot function to order samples
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

DotPlot.2 <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, 
                       dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale.by = "radius", 
                       scale.min = NA, scale.max = NA, order.ids = NULL){
  if (! is.null(assay)){
    assay <- DefaultAssay(object = object)
  }
  #assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  data.features <- FetchData(object = object, vars = features)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)
  }
  else {
    object[[group.by, drop = TRUE]]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]]
    if (length(x = unique(x = splits)) > length(x = cols)) {
      stop("Not enought colors for the number of groups")
    }
    cols <- cols[1:length(x = unique(x = splits))]
    names(x = cols) <- unique(x = splits)
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             data.use <- scale(x = data.use)
                             data.use <- MinMax(data = data.use, min = col.min, 
                                                max = col.max)
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (!is.null(x = split.by)) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = rev(x = features))
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (!is.null(x = split.by)) {
    splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
                                      split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
                         2)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
                     no = "colors")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(order.ids)){
    data.plot$id <- factor(data.plot$id, levels=order.ids)
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + theme_cowplot()
  if (!is.null(x = split.by)) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (is.null(x = split.by)) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}

# run a soft filtering on genes for SCENIC
myFilterGenes <- function(tobj, tassay="RNA", tclusts=NULL, tCountsPerGene=3*0.01*ncol(tobj), tminSamples=ncol(tobj)*0.01){
  # get a subset of cells if requested
  tcells <- colnames(tobj)
  if (! is.null(tclusts)){
    tcells <- names(Idents(tobj))[Idents(tobj) %in% tclusts]
  }
  print(paste(length(tcells),"cells","selected.",sep=" "))
  # The first filter, the total number of reads per gene, is meant to remove genes that are most likely unreliable and provide only noise. 
  # The specific value depends on the data set; for the ones used in the original paper they set the thresholds at, 
  # for example, 3 UMI counts (slightly over the median of the nonzero values) multiplied by 1% of the number of cells in the data set 
  # (e.g., in mouse brain: 3 UMI counts × 30 (1% of cells) = minimum 90 counts per gene). 
  total.UMIs.per.gene <- rowSums(tobj[[tassay]]@counts[,tcells])
  tgenes.filt.1 <- names(total.UMIs.per.gene)[total.UMIs.per.gene >= tCountsPerGene]
  # number of cells in which a gene is detected
  number.detected.cells <- apply(tobj[[tassay]]@counts[,tcells], MARGIN=1, FUN=function(x) { sum(x>0) })
  tgenes.filt.2 <- names(number.detected.cells)[number.detected.cells >= tminSamples]
  # get intersecting genes that pass both filters
  return(intersect(tgenes.filt.1, tgenes.filt.2))
}

# load tumor subtype signature genes
my.load.reference.genes <- function(tobj, tassay='SCT', tcells, tmin.cells=0.05, tref.dir, tgene.sets, tout.file=NULL){
  paper.gene.df <- data.frame(gene=character(), geneset=character())
  # read and clean gene sets
  for (geneset in tgene.sets){
    print(geneset)
    # load reference genes
    tgenes.df <- read.table(file.path(tref.dir, paste(geneset,'mouse','txt',sep='.')), header=T, sep='\t', check.names=F, stringsAsFactors=F)
    print(paste(nrow(tgenes.df), 'genes', 'loaded.'))
    # extract valid genes
    tgenes <- intersect(tgenes.df$MGI.symbol, rownames(tobj[[tassay]]@scale.data))
    print(paste(length(tgenes), 'valid', 'genes', 'in', 'our', 'dataset'))
    # remove genes that express in only small fraction of cells
    tgenes.pop <- apply(tobj[[tassay]]@data[tgenes, tcells], MARGIN=1, FUN=function(x) { sum(x>0) / length(x) })
    tgenes.clean <- names(tgenes.pop)[tgenes.pop > tmin.cells]
    print(paste(length(tgenes.clean), 'genes', 'in', paste0('>',tmin.cells*100,'%'), 'cells', 'are', 'selected'))
    # update gene set and annotation rows
    paper.gene.df <- rbind(paper.gene.df, data.frame(gene=tgenes.clean, geneset=rep(geneset, length(tgenes.clean))))
  }
  
  # check for overlapping genes among gene sets
  # remove overlapping genes!!!
  stat.overlapping.genes <- paper.gene.df %>% group_by(gene) %>% summarise_at('geneset', function(x) {length(x)} ) %>% arrange(desc(geneset))
  overlapping.genes <- subset(stat.overlapping.genes, geneset > 1)[,'gene',drop=T]
  print(paste(length(overlapping.genes), 'overlapping', 'genes', 'removed.'))
  
  paper.gene.df <- paper.gene.df %>% filter(! gene %in% overlapping.genes)
  # write to file
  if (! is.null(tout.file)){
    write.table(paper.gene.df, tout.file, quote=F, sep='\t', row.names=F)
  }
  return(paper.gene.df)
}

# draw heatmap after manually order cells
# version .2 --> order cells (columns) by group
my.manual.heatmap.2 <- function(texp, tanncols, tanncolor, tgroup_column_by, tgroup_column_order=NULL, tcolor, 
                                tquantile=0.1, tclustering_distance_rows="euclidean", tclustering_distance_cols="euclidean", tclustering_method='complete', 
                                tshow_rownames=T, tfig.dir, tsuffix, twidth, theight, tunits, tres, tfontsize_row, tsvg=FALSE){
  # prepare matrix for pheatmap
  texp <- texp[, rownames(tanncols)]
  # set minimum/maximum display value (clip all values that are out of range, 90%/10% quantile as upper/lower bounds)
  tmax <- quantile(c(texp), c(1-tquantile))
  tmin <- quantile(c(texp), c(tquantile))
  texp <- apply(texp, 2, function(x) ifelse(x > tmax, tmax, x))
  texp <- apply(texp, 2, function(x) ifelse(x < tmin, tmin, x))
  # prepare cell (column) order
  cat('calculating column order...')
  if (is.null(tgroup_column_order)){
    tgroup_column_order <- as.vector(unique(tanncols[,tgroup_column_by]))
  }
  tcolumn.order <- c()
  tgaps_col <- c()
  for (tgroup_column in tgroup_column_order){
    tgroup.cells <- tanncols %>% rownames_to_column('cellID') %>% filter(!!sym(tgroup_column_by) %in% tgroup_column) %>% pull('cellID')
    p <- pheatmap(mat=texp[, tgroup.cells], scale="none", cluster_rows=F, cluster_cols=T, clustering_distance_rows=tclustering_distance_rows, 
                  clustering_distance_cols=tclustering_distance_cols, clustering_method=tclustering_method, silent=T)
    tcolumn.order <- c(tcolumn.order, p$tree_col$labels[p$tree_col$order])
    # split columns by group
    if (length(tgaps_col) == 0){
      tgaps_col <- length(tgroup.cells)
    } else {
      tgaps_col <- c(tgaps_col, last(tgaps_col) + length(tgroup.cells))
    }
  }
  cat('ok.\n')
  # draw heatmap plot (no clustering on cells)
  png(file=file.path(tfig.dir,paste('heatmap','enrichment','score',tsuffix,'png',sep='.')), 
      width=twidth, height=theight, units=tunits, pointsize=8, res=tres)
  print(pheatmap(texp[, tcolumn.order], scale="none", color=tcolor, cluster_rows=F, cluster_cols=F, 
                 annotation_col=tanncols, annotation_color=tanncolor, 
                 show_colnames=F, show_rownames=tshow_rownames, gaps_col=tgaps_col, fontsize_row=tfontsize_row))
  dev.off()
  if (tsvg){
    svglite(file=file.path(tfig.dir,paste('heatmap','enrichment','score',tsuffix,'svg',sep='.')), width=twidth, height=theight, pointsize=8)
    print(pheatmap(texp[, tcolumn.order], scale="none", color=tcolor, cluster_rows=F, cluster_cols=F, 
                   annotation_col=tanncols, annotation_color=tanncolor, 
                   show_colnames=F, show_rownames=tshow_rownames, gaps_col=tgaps_col, fontsize_row=tfontsize_row))
    dev.off()
  }
}

# calculate enrichment score on a list of gene sets, and draw heatmap plot to compare enrichment and clusters
my.reference.enrichment.2 <- function(tobj, tassay='SCT', tcells, tmin.cells=0.05, tquantile=0.1, tanncols, tgroup_column_by, tgroup_column_order, 
                                      tclustering_distance_rows="euclidean", tclustering_distance_cols="euclidean", tclustering_method='complete', 
                                      tshow_rownames=T, tref.dir, tgene.sets, tanncolor, tfig.dir, tinfo.dir, tsuffix, tseed=98,
                                      tcolor=colorRampPalette(rev(brewer.pal(n=7, name='RdBu')))(100), 
                                      twidth=15, theight=3.5, tunits="in", tfontsize_row=8, tres=600, tsvg=FALSE){
  # read and clean gene sets
  paper.gene.df <- my.load.reference.genes(tobj=tobj, tassay=tassay, tcells=tcells, tmin.cells=tmin.cells, 
                                           tref.dir=tref.dir, tgene.sets=tgene.sets)
  # create gene set list for enrichment analysis
  paper.gene.list <- list()
  for (gset in tgene.sets){
    tdf <- paper.gene.df %>% filter(geneset == gset)
    paper.gene.list[[gset]] <- tdf[, 'gene', drop=T]
  }
  # number of genes per set
  print(sapply(paper.gene.list, FUN=length))
  # run AddModuleScore to score each gene set
  tenrich.name <- 'GSE'
  tctrl = min(vapply(X=paper.gene.list, FUN=length, FUN.VALUE=numeric(length=1)))
  #print(tctrl)
  tobj.tmp <- AddModuleScore(object=subset(tobj, cells=tcells), features=paper.gene.list, pool=NULL, name=tenrich.name, ctrl=tctrl, 
                             assay=tassay, seed=tseed)
  # collect enrichment score for each gene set
  tscores <- FetchData(tobj.tmp, vars=c(grep(tenrich.name, colnames(tobj.tmp@meta.data), value=T)))
  colnames(tscores) <- names(paper.gene.list)
  # write per-cell enrichment scores to file
  write.table(tscores, file=file.path(tinfo.dir,paste("enrichment","scores",tsuffix,"txt",sep=".")), quote=F, sep="\t", row.names=T, col.names=NA)
  # draw heatmap plot
  my.manual.heatmap.2(texp=t(tscores), tanncols=tanncols, tanncolor=tanncolor, 
                      tgroup_column_by=tgroup_column_by, tgroup_column_order=tgroup_column_order, tcolor=tcolor, tquantile=tquantile, 
                      tclustering_distance_rows=tclustering_distance_rows, tclustering_distance_cols=tclustering_distance_cols, tclustering_method=tclustering_method, 
                      tshow_rownames=tshow_rownames, tfig.dir=tfig.dir, tsuffix=tsuffix, twidth=twidth, theight=theight, 
                      tunits=tunits, tres=tres, tfontsize_row=tfontsize_row, tsvg=tsvg)
  # draw violin plot
  tdata <- tscores %>% rownames_to_column('cellID') %>% left_join(tanncols %>% rownames_to_column('cellID'), by='cellID') %>% 
    gather(key="subtype", value="score", tgene.sets)
  g <- ggplot(tdata, aes_string(x=tgroup_column_by, y='score', fill=tgroup_column_by)) + geom_violin(scale="width")
  g <- g + scale_fill_manual(values=tanncolor[[tgroup_column_by]])
  g <- g + stat_summary(fun.data=median_hilow, geom="pointrange", color="gray85")
  g <- g + geom_hline(yintercept=0, linetype="dashed", color = "red")
  g <- g + coord_flip() + facet_wrap(~subtype) + theme_bw()
  ggsave(file.path(tfig.dir, paste('violin','enrichment','scores',tsuffix,'png', sep='.')), width=6.5, height=6.5, dpi=300)
  if (tsvg){
    ggsave(file.path(tfig.dir, paste('violin','enrichment','scores',tsuffix,'svg', sep='.')), width=6.5, height=6.5)
  }
}

# Heatmap plot showing gene expressions in cells
# reorder rows and columns by groups such that items within the same group stay together
# !!! cells extracted based on tanncols data.frame
# !!! genes extracted based on trowcols data.frame
myHeatmap.3 <- function(tobj, toutfile, tassay="RNA", 
                        tanncols, tannrows, tanncolor=NA, 
                        tgroup_row_by, tgroup_row_order=NULL, tgroup_column_by, tgroup_column_order=NULL,
                        tclustering_distance_rows="euclidean", tclustering_distance_cols="euclidean", tclustering_method='complete', tshow_rownames=T, 
                        tmax=3, tmin=-3, tcolor=colorRampPalette(rev(brewer.pal(n=7, name='RdBu')))(100), 
                        twidth=8, theight=10, tunits="in", tfontsize_row=8, tres=300, tsvg=FALSE){
  # valid genes?
  tgenes <- intersect(rownames(tannrows), rownames(tobj))
  if (length(tgenes) > 0){
    cat(paste(length(tgenes),'genes to plot.\n',sep=' '))
  } else {
    cat('No valid genes found.\n')
    return(NULL)
  }
  # get cells
  tcells <- rownames(tanncols)
  # prepare matrix for pheatmap
  texp <- tobj[[tassay]]@scale.data[tgenes, tcells]
  # set minimum/maximum display value (clip all values that are out of range)
  texp <- apply(texp, 2, function(x) ifelse(x > tmax, tmax, x))
  texp <- apply(texp, 2, function(x) ifelse(x < tmin, tmin, x))
  # prepare gene (row) order
  cat('calculating row order...')
  if (is.null(tgroup_row_order)){
    tgroup_row_order <- as.vector(unique(tannrows[,tgroup_row_by]))
  }
  trow.order <- c()
  tgaps_row <- c()
  for (tgroup_row in tgroup_row_order){
    tgroup.genes <- tannrows %>% rownames_to_column('gene') %>% filter(!!sym(tgroup_row_by) %in% tgroup_row) %>% pull('gene')
    p <- pheatmap(mat=texp[tgroup.genes, ], scale="none", cluster_rows=T, cluster_cols=T, clustering_distance_rows=tclustering_distance_rows, 
                  clustering_distance_cols=tclustering_distance_cols, clustering_method=tclustering_method, silent=T)
    trow.order <- c(trow.order, p$tree_row$labels[p$tree_row$order])
    # split rows by group
    if (length(tgaps_row) == 0){
      tgaps_row <- length(tgroup.genes)
    } else {
      tgaps_row <- c(tgaps_row, last(tgaps_row) + length(tgroup.genes))
    }
  }
  cat('ok.\n')
  # prepare cell (column) order
  cat('calculating column order...')
  if (is.null(tgroup_column_order)){
    tgroup_column_order <- as.vector(unique(tanncols[,tgroup_column_by]))
  }
  tcolumn.order <- c()
  tgaps_col <- c()
  for (tgroup_column in tgroup_column_order){
    tgroup.cells <- tanncols %>% rownames_to_column('cellID') %>% filter(!!sym(tgroup_column_by) %in% tgroup_column) %>% pull('cellID')
    p <- pheatmap(mat=texp[trow.order, tgroup.cells], scale="none", cluster_rows=F, cluster_cols=T, clustering_distance_rows=tclustering_distance_rows, 
                  clustering_distance_cols=tclustering_distance_cols, clustering_method=tclustering_method, silent=T)
    tcolumn.order <- c(tcolumn.order, p$tree_col$labels[p$tree_col$order])
    # split columns by group
    if (length(tgaps_col) == 0){
      tgaps_col <- length(tgroup.cells)
    } else {
      tgaps_col <- c(tgaps_col, last(tgaps_col) + length(tgroup.cells))
    }
  }
  cat('ok.\n')
  # plot
  png(file=toutfile, width=twidth, height=theight, units=tunits, pointsize=8, res=tres)
  print(pheatmap(texp[trow.order, tcolumn.order], scale="none", color=tcolor, cluster_rows=F, cluster_cols=F, 
                 annotation_col=tanncols, annotation_row=tannrows, annotation_color=tanncolor, 
                 show_colnames=F, show_rownames=tshow_rownames, gaps_col=tgaps_col, gaps_row=tgaps_row, fontsize_row=tfontsize_row))
  dev.off()
  if (tsvg){
    svglite(file=gsub(file_ext(toutfile), 'svg', toutfile), width=twidth, height=theight, pointsize=8)
    print(pheatmap(texp[trow.order, tcolumn.order], scale="none", color=tcolor, cluster_rows=F, cluster_cols=F, 
                   annotation_col=tanncols, annotation_row=tannrows, annotation_color=tanncolor, 
                   show_colnames=F, show_rownames=tshow_rownames, gaps_col=tgaps_col, gaps_row=tgaps_row, fontsize_row=tfontsize_row))
    dev.off()
  }
}
