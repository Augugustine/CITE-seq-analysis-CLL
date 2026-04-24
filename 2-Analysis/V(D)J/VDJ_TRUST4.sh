#!/bin/bash
#SBATCH -J script
#SBATCH -o VDJ_trust4.out
#SBATCH -e VDJ_trust4_err.out
#SBATCH -t 24:00:00
#SBATCH --mem=16G
#SBATCH -c 2

 
module load devel/Miniconda/Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
conda activate trust4_env

for i in 1 4 8 11 14; do
  run-trust4 \
  -t 2 \
  -b /home/ablanc/work/NetBIO2/CNA/Numbat_p2/p2_D${i}.bam \
  -f /work/user/ablanc/NetBIO2/VDJ/hg38_bcrtcr.fa \
  --ref /work/user/ablanc/NetBIO2/VDJ/human_IMGT+C.fa \
  --barcode CB \
  --UMI UB \
  --od vdj_trust4_D${i} \
  -o vdj_p2_D${i}
  
done


