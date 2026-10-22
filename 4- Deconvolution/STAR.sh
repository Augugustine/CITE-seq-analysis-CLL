#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --error=logs/star_%A_%a.err
#SBATCH --array=0-57                   
#SBATCH --cpus-per-task=12         
#SBATCH --mem=40G
#SBATCH --time=48:00:00 

# Module
module load bioinfo/STAR/2.7.11b
module load bioinfo/fastp/1.3.6
#module load zlib

# Folders
GENOME_DIR="STAR_index/STAR_index"
FASTQ_DIR="FastQ"
OUT_DIR="Count_results"

mkdir -p $OUT_DIR
mkdir -p logs

# List samples
cd $FASTQ_DIR
samples=( $(ls *_1.fastq.gz | sed 's/_1\.fastq\.gz//') )
cd ..

# Array task
current_sample=${samples[$SLURM_ARRAY_TASK_ID]}

echo "Fastp: $current_sample"

echo "Trimming: $current_sample"
fastp -i ${FASTQ_DIR}/${current_sample}_1.fastq.gz -I ${FASTQ_DIR}/${current_sample}_2.fastq.gz \
      -o ${FASTQ_DIR}/${current_sample}_clean_R1.fastq.gz -O ${FASTQ_DIR}/${current_sample}_clean_R2.fastq.gz \
      --thread 6

echo "STAR: $current_sample"

# STAR
STAR --runThreadN 12 \
     --genomeDir $GENOME_DIR \
     --readFilesIn ${FASTQ_DIR}/${current_sample}_clean_R1.fastq.gz ${FASTQ_DIR}/${current_sample}_clean_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUT_DIR}/${current_sample}_ \
     --outSAMtype None \
     --quantMode GeneCounts

#rm ${FASTQ_DIR}/${current_sample}_clean_R1.fastq.gz ${FASTQ_DIR}/${current_sample}_clean_R2.fastq.gz

echo "End for sample $current_sample"
