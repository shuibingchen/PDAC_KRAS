#!/usr/bin/env python
# process_scenic_results.py
# further process results from pySCENIC CLI commands
# 

import os.path
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
import pandas as pd

workdir = "."
infodir = os.path.join(workdir, "info")

aucMtxFile = os.path.join(infodir, "auc_mtx.csv")
binMtxFile = os.path.join(infodir, "bin_mtx.csv")
thrFile = os.path.join(infodir, "threshold.csv")
umapFile = os.path.join(infodir, "umap.auc_mtx.txt")

# Regulon activity binarization
auc_mtx = pd.read_csv(aucMtxFile, index_col=0)
#print(auc_mtx.shape)
# the input auc matrix should be nCells x nRegulons
bin_mtx, thresholds = binarize(auc_mtx.T)
#print(bin_mtx.shape)
#print(thresholds.shape)

# write to file
bin_mtx.to_csv(binMtxFile)
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(thrFile)

print("Complete!")

