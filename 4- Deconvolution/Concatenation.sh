#!/bin/bash
#SBATCH --job-name=concatenation
#SBATCH --output=concatenation.out
#SBATCH --error=concatenation.err            
#SBATCH --mem=8G
#SBATCH --time=04:00:00 

files=(*_ReadsPerGene.out.tab)

# Header
{
  printf "GeneID"
  for f in "${files[@]}"; do printf "\t%s" "${f%_ReadsPerGene.out.tab}"; done
  printf "\n"
} > count_matrix.txt

# Body: one row per gene, one column per sample (same order as the header)
awk 'BEGIN{FS=OFS="\t"}
     FNR>4 {                      # skip the 4 N_ summary lines
       if (FNR==NR) id[FNR]=$1    # gene IDs taken from the first file
       row[FNR] = row[FNR] OFS $4 # append column 4 of each file
       n = FNR
     }
     END { for (i=5; i<=n; i++) print id[i] row[i] }' "${files[@]}" >> count_matrix.txt

 
# Verification 

head -c 300 count_matrix.txt | head -3
awk -F'\t' '{print NF}' count_matrix.txt | sort -u   # nb samples + 1
wc -l count_matrix.txt                                # = nb genes + 1

# Controle in 1 sample
sample=$(ls *_ReadsPerGene.out.tab | sed -n 2p)
diff <(tail -n +5 "$sample" | cut -f4) <(tail -n +2 count_matrix.txt | cut -f3) && echo OK
