# step_5.generate_plots.R
# step 5: generate plots in the manuscript
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

# set color for clusters
my.cluster.color <- c(brewer.pal(6,'Paired'))
names(my.cluster.color) <- 0:5

# set color for samples
my.sample.color <- brewer.pal(10, 'Paired')
names(my.sample.color) <- c("WT","WT-P","KC","KC-P","KPC","KPC-P","KSC","KSC-P","KPSC","KPSC-P")

# set color for gene expressions
myExpLowColor <- '#d9d9d9'
myExpHighColor <- '#b30000'

# set x/y-axis boundaries of UMAP plot
x.upper <- 14
x.lower <- -11
y.upper <- 13
y.lower <- -8

# retrieve cells in untreated samples
untreated.samples <- grep('-P', sample.info$Name, invert=T, value=T)
untreated.cells <- rownames(FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') 
                            %>% filter(Name %in% untreated.samples) %>% column_to_rownames('cellID'))

# Figure 1B: UMAP plot illustrating cell clusters in untreated samples
g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", tcells=untreated.cells, 
               tcolor=my.cluster.color, tlabel=TRUE, tsplit=FALSE, txlim=c(x.lower,x.upper), 
               tylim=c(y.lower,y.upper), tptsize=1.5, talpha=0.7, tltsize=20, tatlsize=22)
ggsave(file.path(figdir, "Figure.1B.png"), plot=g, width=10, height=8, dpi=300)

# Figure 1C: UMAP plot illustrating sample classifications based on CMO-labels
g <- myDimPlot(tobj=panc, treduct="umap", tcate="Name", tsuffix="Sample", tcells=untreated.cells, 
               tcolor=my.sample.color, tlabel=FALSE, tsplit=FALSE, txlim=c(x.lower,x.upper), 
               tylim=c(y.lower,y.upper), tptsize=1.5, talpha=0.7, tltsize=20, tatlsize=22)
ggsave(file.path(figdir, "Figure.1C.png"), plot=g, width=10, height=8, dpi=300)

# Figure S7A: UMAP plot illustrating cell clusters in untreated and treated samples
g <- myDimPlot3(tobj=panc, treduct="umap", tgroup_by="ident", tgroup_order=0:5, 
                tsuffix="Treatment", tcolor=my.cluster.color, tsplit_by='Treatment',
                tsplit_order=c('Untreated','Treated'), tlabel=T, tncol=2, tptsize=0.5, tlbsize=4) + theme(legend.position="none")
ggsave(file.path(figdir, "Figure.S7A.png"), plot=g, width=9, height=5, dpi=300)

# Figure 1E: UMAP plot showing the expression of pancreatic progenitor pattern genes, including Sox9, Foxa2, Hnf1β, Gata6 and Pdx1
g <- MyFeaturePlot(tobj=panc, tgenes=c('Sox9','Foxa2','Hnf1b','Gata6','Pdx1'), tcells=untreated.cells, tassay="SCT", 
                   treduction.name="umap", tncol=5, tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tlegend=NULL)
ggsave(file.path(figdir, "Figure.1E.png"), plot=g, width=21, height=3.5, dpi=300)

# Figure S1F: UMAP plot showing the expression of acinar cell (Gata4 and Ptf1a) and endocrine lineage (Nkx2.2, Nkx6.1 and Neurog3) marker genes.
g <- MyFeaturePlot(tobj=panc, tgenes=c('Gata4','Nkx6-1'), tcells=untreated.cells, tassay="SCT", 
                   treduction.name="umap", tncol=2, tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tlegend=NULL)
ggsave(file.path(figdir, "Figure.S1F.a.png"), plot=g, width=8, height=3.5, dpi=300)
# no cells expressing Neurog3
# use 'RNA' slot to plot this gene since 'SCT' automatically remove genes with very few cells expressing it
g <- MyFeaturePlot(tobj=panc, tgenes='Neurog3', tcells=untreated.cells, tassay="RNA", 
                   treduction.name="umap", tncol=1, tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tlegend=NULL)
ggsave(file.path(figdir, "Figure.S1F.b.png"), plot=g, width=4, height=3.5, dpi=300)
# very few cells expressing Nkx2-2 and Ptf1a
# use 'RNA' slot to plot this gene since 'SCT' automatically remove genes with very few cells expressing it
g <- MyFeaturePlot(tobj=panc, tgenes=c('Nkx2-2','Ptf1a'), tcells=untreated.cells, tassay="RNA", 
                   treduction.name="umap", tncol=1, tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tlegend=NULL)
ggsave(file.path(figdir, "Figure.S1F.c.png"), plot=g, width=8, height=3.5, dpi=300)

# Figure 1D: Correlation analysis of cell clusters based on unsupervised analysis and CMO-labels.
corr.clusters.samples <- as.matrix(as.data.frame.matrix(table(FetchData(panc, vars=c('Name','ident')))) %>%
                                     rownames_to_column('sample') %>% gather(key='cluster',value='ncells',as.character(0:5)) %>%
                                     group_by(sample) %>% mutate(pcells=ncells/sum(ncells)) %>% select(sample, cluster, pcells) %>%
                                     spread(cluster, pcells) %>% column_to_rownames('sample'))

untreated.samples.order <- c('WT','KC','KPC','KSC','KPSC')

png(file.path(figdir, 'Figure.1D.png'), width=6, height=2.8, units='in', res=300)
pheatmap(corr.clusters.samples[untreated.samples.order,], color=colorRampPalette(rev(brewer.pal(n=7, name='RdBu')))(100),
         scale="none", cluster_cols=F, cluster_rows=F, angle_col=0, fontsize=16,
         breaks=seq(0,1,0.01), legend_breaks=c(0,0.2,0.4,0.6,0.8,1), legend_labels=as.character(c(0,0.2,0.4,0.6,0.8,1)))
dev.off()

# Figure 1G: Dot plot comparing the relative expression levels for Gata6 of KPC, KSC and KPSC organoids.
untreated.KPC.KSC.KPSC.cells <- rownames(FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>% 
                                           filter(Name %in% c('KPC','KSC','KPSC')) %>% column_to_rownames('cellID'))
g <- DotPlot.2(subset(panc, cells=untreated.KPC.KSC.KPSC.cells), assay='SCT', features='Gata6', group.by='Name', 
             order.ids=c('KPC','KSC','KPSC'))
ggsave(file.path(figdir, 'Figure.1G.png'), plot=g, width=5.5, height=4, dpi=300)

# Figure 5C: Dot plots describing the relative expression levels for genes of the cholesterol biosynthesis pathway in WT, KC, KPC, KSC, and KPSC organoids
cholesterol.biosynthesis.genes <- c("Acat2","Apoe","Lbr","Hmgcs1","Hmgcr","Lss","Mvd",
                                    "Pmvk","Msmo1","Hsd17b7","Dhcr7","Sc5d","Ldlr","Vldlr")

g <- DotPlot.2(subset(panc, cells=untreated.cells), assay='SCT', features=cholesterol.biosynthesis.genes, 
               cols="RdBu", group.by='Name', order.ids=c('WT','KC','KPC','KSC','KPSC')) + 
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Figure.5C.png"), plot=g, width=4.5, height=6, dpi=300)

# Figure 5I-L: Dot plots comparing the relative expression levels for genes of the cholesterol biosynthesis pathway of 
#              Perhexiline maleate or DMSO treated KC (I), KPC (J), KSC (K) and KPSC (L) organoids. 
#              P: Perhexiline maleate; C: control.
KC.cells <- rownames(FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>%
                       filter(Name %in% c('KC','KC-P')) %>% column_to_rownames('cellID'))
KPC.cells <- rownames(FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>%
                        filter(Name %in% c('KPC','KPC-P')) %>% column_to_rownames('cellID'))
KSC.cells <- rownames(FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>%
                        filter(Name %in% c('KSC','KSC-P')) %>% column_to_rownames('cellID'))
KPSC.cells <- rownames(FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>%
                         filter(Name %in% c('KPSC','KPSC-P')) %>% column_to_rownames('cellID'))

g <- DotPlot(subset(panc, cells=KC.cells), assay='SCT', features=cholesterol.biosynthesis.genes, cols="RdBu", group.by='Name') +
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Figure.5I.png"), plot=g, width=3.5, height=6.5, dpi=300)

g <- DotPlot(subset(panc, cells=KPC.cells), assay='SCT', features=cholesterol.biosynthesis.genes, cols="RdBu", group.by='Name') +
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Figure.5J.png"), plot=g, width=3.5, height=6.5, dpi=300)

g <- DotPlot(subset(panc, cells=KSC.cells), assay='SCT', features=cholesterol.biosynthesis.genes, cols="RdBu", group.by='Name') +
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Figure.5K.png"), plot=g, width=3.5, height=6.5, dpi=300)

g <- DotPlot(subset(panc, cells=KPSC.cells), assay='SCT', features=cholesterol.biosynthesis.genes, cols="RdBu", group.by='Name') +
  theme_bw() + coord_flip() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggsave(file.path(figdir, "Figure.5L.png"), plot=plot, width=3.5, height=6.5, dpi=300)


# get cells from KPC, KSC and KPSC samples
untreated.KPC.KSC.KPSC.cells <- rownames(FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>% 
                                           filter(Name %in% c('KPC','KSC','KPSC')) %>% column_to_rownames('cellID'))

# set color for sample KPC, KSC and KPSC
my.KPC.KSC.KPSC.color <- c(brewer.pal(9,'Set1'))[2:4]
names(my.KPC.KSC.KPSC.color) <- c('KPC','KSC','KPSC')

# load binarized regulon activity
binReg <- read.table(file.path(infodir, "bin_mtx.csv"), header=T, check.names=F, stringsAsFactors=F, sep=",", row.names=1)

# prepare cell annotations in heatmap
anncols <- FetchData(panc, vars=c("Name"), cells=untreated.KPC.KSC.KPSC.cells)

# set color scheme in heatmap plot
my.color <- colorRampPalette(c("grey95", "grey50", "grey5"))(100)

# Figure 1F: Heatmap showing the master regulators of KSC, KPC and KPSC organoids
png(file.path(figdir, "Figure.1F.pngg"), width=16, height=12.5, units="in", res=300)
pheatmap(as.matrix(t(binReg[untreated.KPC.KSC.KPSC.cells,])), 
         annotation_col=anncols, color=my.color, scale="none",
         annotation_colors=list(Name=my.KPC.KSC.KPSC.color),
         cluster_rows=T, cluster_cols=F, show_colnames=F, show_rownames=T, fontsize_col=12, fontsize_row=3)
dev.off()

# Figure S2C-D: Tumor subtypes (basal and classical) scores of KPC, KSC, and KPSC organoids. 
untreated.KPC.KSC.KPSC.cells <- rownames(FetchData(panc, vars=c('Name')) %>% rownames_to_column('cellID') %>% 
                                           filter(Name %in% c('KPC','KSC','KPSC')) %>% column_to_rownames('cellID'))

# load basal and classical signature genes
sig.genes.ref_1 <- my.load.reference.genes(tobj=panc, tassay='SCT', tcells=untreated.KPC.KSC.KPSC.cells, tmin.cells=0.01,
                                           tref.dir=file.path(sourcedir,'subtype','ref_1'),
                                           tgene.sets=c('Basal','Classical'),
                                           tout.file=file.path(sourcedir,'subtype','ref_1','cleaned.gene.set.txt'))

# prepare gene annotations
annrows.manual.ref_1 <- sig.genes.ref_1 %>% column_to_rownames('gene')

# prepare cell annotations in heatmap plot
anncols.manual <- FetchData(panc, cells=untreated.KPC.KSC.KPSC.cells, vars=c('Name'))
anncols.manual$Name <- factor(anncols.manual$Name, levels=c('KPC','KSC','KPSC'))

# calculate enrichment scores
my.reference.enrichment.2(tobj=panc, tassay='SCT', tcells=untreated.KPC.KSC.KPSC.cells, tmin.cells=0.01, tquantile=0.1, 
                          tanncols=anncols.manual, tgroup_column_by='Name', tgroup_column_order=c('KPC','KSC','KPSC'),
                          tclustering_distance_rows="euclidean", tclustering_distance_cols="euclidean", tclustering_method='complete',
                          tshow_rownames=T, tref.dir=file.path(sourcedir,'subtype','ref_1'), tgene.sets=c('Classical','Basal'),
                          tanncolor=list(Name=my.KPC.KSC.KPSC.color), tfig.dir=figdir, tinfo.dir=infodir,
                          tsuffix='Figure.S2C-D', tseed=98,
                          tcolor=colorRampPalette(rev(brewer.pal(n=7, name='RdBu')))(100),
                          twidth=15, theight=2.5, tunits="in", tfontsize_row=8, tres=600, tsvg=F)

# Figure S2E: Heatmap showing the expression of Moffitt subtype signature genes across individual clusters of KPC, KSC, and KPSC organoids.
# order cells and columns by group
myHeatmap.3(tobj=panc,
            toutfile=file.path(figdir,'Figure.2E.png'),
            tassay="SCT", tanncols=anncols.manual, tannrows=annrows.manual.ref_1,
            tanncolor=list(Name=my.KPC.KSC.KPSC.color, geneset=c('Basal'='#377eb8','Classical'='#4daf4a')),
            tgroup_row_by='geneset', tgroup_row_order=c('Classical','Basal'),
            tgroup_column_by='Name', tgroup_column_order=c('KPC','KSC','KPSC'),
            tclustering_distance_rows="euclidean", tclustering_distance_cols="euclidean", tclustering_method='complete', tshow_rownames=T,
            tmax=3, tmin=-3, tcolor=colorRampPalette(rev(brewer.pal(n=7, name='RdBu')))(100),
            twidth=10.5, theight=7.5, tunits="in", tfontsize_row=12, tres=300, tsvg=F)

# Figure S2F-G: Tumor subtypes (sig.1, sig.2, sig.6 and sig.10) scores of KPC, KSC, and KPSC organoids.

# load basal-A/B and classical-A/B signature genes
sig.genes.ref_2 <- my.load.reference.genes(tobj=panc, tassay='SCT', tcells=untreated.KPC.KSC.KPSC.cells, tmin.cells=0.01, 
                                           tref.dir=file.path(sourcedir,'subtype','ref_2'), 
                                           tgene.sets=c('sig.1','sig.2','sig.6','sig.10'), 
                                           tout.file=file.path(sourcedir,'subtype','ref_2','cleaned.gene.set.ref_2.txt'))

# calculate enrichment scores
my.reference.enrichment.2(tobj=panc, tassay='SCT', tcells=untreated.KPC.KSC.KPSC.cells, tmin.cells=0.01, tquantile=0.1, 
                          tanncols=anncols.manual, tgroup_column_by='Name', tgroup_column_order=c('KPC','KSC','KPSC'), 
                          tclustering_distance_rows="euclidean", tclustering_distance_cols="euclidean", tclustering_method='complete', 
                          tshow_rownames=T, tref.dir=file.path(sourcedir,'subtype','ref_2'), 
                          tgene.sets=c('sig.1','sig.2','sig.6','sig.10'), 
                          tanncolor=list(Name=my.KPC.KSC.KPSC.color), tfig.dir=figdir, tinfo.dir=infodir, 
                          tsuffix='Figure.S2F-G', tseed=98,
                          tcolor=colorRampPalette(rev(brewer.pal(n=7, name='RdBu')))(100), 
                          twidth=15, theight=4, tunits="in", tfontsize_row=8, tres=600, tsvg=F)
