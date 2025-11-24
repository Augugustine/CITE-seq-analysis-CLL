#!/bin/bash
#SBATCH -o pyscenic_hsc_sev.out
#SBATCH -e pyscenic_hsc_sev.err
#SBATCH -J scenic_HSC_SEV
#SBATCH --nodes=1               
#SBATCH -n 20
#SBATCH --mem=200G
#SBATCH -t 48:00:00

# Initialiser conda
source $(conda info --base)/etc/profile.d/conda.sh

# Activer votre environnement
conda activate scenic_env

# human
f_db_names="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
f_motif_path="motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
f_tf_list="allTFs_hg38.txt"

dir_result="results/"
input_loom="TimeCourse.loom"


# Step 1
echo "Step 1 pyscenic grn start - $(date)"
arboreto_with_multiprocessing.py \
    ${input_loom} \
    ${f_tf_list} \
    --method grnboost2 \
    --output ${dir_result}/step_1_fibo_grn.tsv \
    --num_workers 20 \
    --seed 777
echo "Step 1 pyscenic grn finish - $(date)"

# Step 2
echo "Step 2 pyscenic ctx start - $(date)"
pyscenic ctx ${dir_result}/step_1_fibo_grn.tsv \
 ${f_db_names} \
--annotations_fname ${f_motif_path} \
--expression_mtx_fname ${input_loom} \
--output ${dir_result}/step_2_reg.csv \
--mask_dropouts \
--num_workers 20
echo "Step 2 pyscenic ctx finish - $(date)"

# Step 3
echo "Step 3 pyscenic aucell start - $(date)"
pyscenic aucell \
 ${input_loom} \
 ${dir_result}/step_2_reg.csv \
--seed 21 \
--output ${dir_result}/step_3_aucell.csv \
--num_workers 20
echo "Step 3 pyscenic aucell finish - $(date)"

echo "All finish - $(date)"
