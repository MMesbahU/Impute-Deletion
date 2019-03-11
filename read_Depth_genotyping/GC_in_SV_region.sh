#!/bin/bash

# input: list of deletions (or SVs) in tab separated. 1st 3 columns should be: Chr del_Start del_END
# The fasta index file should be in the same directory
# to make index (uncompressed or compressed by bgzip): samtools faidx ${Reference}
Reference=${1}
sv_list=${2}
#
bedtools nuc -fi ${Reference} -bed ${sv_list} > GC_content_in_SV_regions.txt

# Output file: tab sepatered file with header
## 12 columns in "GC_content_in_SV_regions.txt"
## Chrom[1]:Del_Start[2]:Del_End[3]:AT_pct[4]:GC_pct[5]:num_A[6]:num_C[7]:num_G[8]:num_T[9]:num_N[10]:num_oth[11]:seq_len[12] 
## GC-conten in column 5

awk 'NR>1 {printf "%s\t%d\t%d\t%0.0f\n", $1, $2, $3, ($5*100)} ' GC_content_in_SV_regions

   

