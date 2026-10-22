#!/bin/bash
#SBATCH --job-name=download
#SBATCH --output=download.out
#SBATCH --error=download.err               
#SBATCH --cpus-per-task=4       
#SBATCH --mem=64G
#SBATCH --time=24:00:00 

# Module
module load bioinfo/SRA-Toolkit/3.3.0

#prefetch --option-file SRA.txt

mkdir -p FastQ_final
mkdir -p logs

SRA_FILE="SRA.txt"

while read -r SRR; do
    [[ -z "$SRR" || "$SRR" =~ ^# ]] && continue

    echo "Sample: $SRR"
    sra_path="FastQ/${SRR}"

    #fasterq-dump
    fasterq-dump "$sra_path" \
        --split-files \
        --threads 4 \
        --outdir FastQ_final

    # Compression
    echo "Compression FASTQ sample $SRR..."
    gzip -f FastQ_final/${SRR}_1.fastq
    gzip -f FastQ_final/${SRR}_2.fastq

done < "$SRA_FILE"
