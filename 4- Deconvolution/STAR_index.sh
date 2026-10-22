#!/bin/bash
#SBATCH -J star
#SBATCH -o output.out
#SBATCH -e error.out
#SBATCH -t 04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12


# Module
module load bioinfo/STAR/2.7.11b

# download fasta and gtf files, GRCh38
#wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.116.gtf.gz


STAR --runMode genomeGenerate \
    --runThreadN 12 \
    --genomeDir STAR_index \
    --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile Homo_sapiens.GRCh38.116.gtf \
    --sjdbOverhang 99


