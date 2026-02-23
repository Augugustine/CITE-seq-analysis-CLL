#!/bin/bash
#SBATCH -J script
#SBATCH -o output.out
#SBATCH -e error.out
#SBATCH -t 48:00:00
#SBATCH --mem=124G
#SBATCH -c 16

# Les barcodes dans les fichiers ne doivent pas etre entourés de ""
#for i in 4 8 11 14; do
#    sed 's/"//g' Numbat/barcodes_D${i}.tsv > Numbat/barcodes_D${i}_clean.tsv
#    rm Numbat/barcodes_D${i}.tsv
#    mv Numbat/barcodes_D${i}_clean.tsv Numbat/barcodes_D${i}.tsv
#done


module load containers/singularity/3.9.9

WORKDIR=/home/ablanc/work/NetBIO2/CNA/Numbat
IMG=$WORKDIR/numbat-rbase_latest.sif

singularity exec \
  -B $WORKDIR:/mnt/mydata \
  -B /home/ablanc/work/NetBIO2/CNA/Numbat \
  $IMG \
  Rscript /numbat/inst/bin/pileup_and_phase.R \
      --label Patient2 \
      --samples D1,D4,D8,D11,D14 \
      --bams /mnt/mydata/p2_D1.bam,/mnt/mydata/p2_D4.bam,/mnt/mydata/p2_D8.bam,/mnt/mydata/p2_D11.bam,/mnt/mydata/p2_D14.bam \
      --barcodes /mnt/mydata/barcodes_D1.tsv,/mnt/mydata/barcodes_D4.tsv,/mnt/mydata/barcodes_D8.tsv,/mnt/mydata/barcodes_D11.tsv,/mnt/mydata/barcodes_D14.tsv \
      --outdir /mnt/mydata/P2_Numbat \
      --gmap /Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
      --snpvcf /data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
      --paneldir /data/1000G_hg38 \
      --ncores 16