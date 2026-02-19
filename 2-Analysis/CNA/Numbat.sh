#!/bin/bash
#SBATCH -J script
#SBATCH -o output.out
#SBATCH -e error.out
#SBATCH -t 48:00:00
#SBATCH --mem=124G
#SBATCH -c 16


module load containers/singularity/3.9.9

WORKDIR=/home/ablanc/work/NetBIO2/CNA/Numbat
IMG=$WORKDIR/numbat-rbase_latest.sif

singularity exec \
  -B $WORKDIR:/mnt/mydata \
  -B /home/ablanc/work/NetBIO2/CNA/Numbat \
  $IMG \
  Rscript /numbat/inst/bin/pileup_and_phase.R \
      --label Patient2 \
      --samples D1\
      --bams /mnt/mydata/p2_D1.bam \
      --barcodes /mnt/mydata/barcodes_D1.tsv \
      --outdir /mnt/mydata/Test \
      --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
      --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
      --paneldir /data/1000G_hg38 \
      --ncores 16