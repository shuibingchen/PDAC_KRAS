#!/usr/bin/env sh
# step_4.reconstruct_gene_regulatory_networks.sh
# step 4: reconstruct gene regulatory networks, create regulons, and quantify regulon activities at a cellular resolution
#

workdir=.
srcdir=${workdir}/src
infodir=${workdir}/info
logdir=${infodir}/logs
dbdir=${workdir}/data/cisTarget

# random seed
rseed=98

# number of cores
ncpu=12

expfile=${infodir}/UMI.counts.cleaned.csv

# use pyscenic v0.12.1

# step 1: network inference
echo -n "step 1: pyscenic grn..."
pyscenic grn \
        --output ${infodir}/expr_mat.adjacencies.tsv \
        --transpose \
        --method grnboost2 \
        --seed ${rseed} \
        --num_workers ${ncpu} \
        ${expfile} \
        ${dbdir}/allTFs_mm.txt \
>${logdir}/pyscenic.grn.log 2>&1
echo "ok."

# step 2: Candidate regulon generation and regulon prediction
echo -n "step 2: pyscenic ctx..."
pyscenic ctx \
        --output ${infodir}/regulons.csv \
        --transpose \
        --annotations_fname ${dbdir}/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
        --num_workers ${ncpu} \
        --expression_mtx_fname ${expfile} \
        --mask_dropouts \
        ${infodir}/expr_mat.adjacencies.tsv \
        ${dbdir}/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
>${logdir}/pyscenic.ctx.log 2>&1
echo "ok."

# step 3: Cellular enrichment 
echo -n "step 3: pyscenic aucell..."
pyscenic aucell \
        --output ${infodir}/auc_mtx.csv \
        --transpose \
        --num_workers ${ncpu} \
        --seed ${rseed} \
        ${expfile} \
        ${infodir}/regulons.csv \
>${logdir}/pyscenic.aucell.log 2>&1
echo "ok."

# step 4: Binarize AUC matrix + Visualize AUC results with UMAP dimensionality reduction
python ${srcdir}/process_scenic_results.py >${logdir}/process_scenic_results.log 2>&1

echo "All Complete!"

