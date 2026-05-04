#!/bin/bash
#SBATCH -J mixcr_fastq
#SBATCH -o VDJ_mixcr.out
#SBATCH -e VDJ_mixcr_err.out
#SBATCH -t 24:00:00
#SBATCH --mem=16G
#SBATCH -c 4

 
module load devel/Miniconda/Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mixcr_env
#mixcr activate-license 

mixcr analyze 10x-sc-5gex \
    --species hsa \
    /home/ablanc/work/NetBIO2/Data/Fastq/P2/Undetermined_S0_L{{n}}_{{R}}_001.fastq.gz \
    p2_mixcr_fastq \
    -t 4 
    


