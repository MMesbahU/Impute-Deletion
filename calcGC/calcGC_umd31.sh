#!/bin/bash

## GC for Bovine Genome Assembly UMD3.1

# To Run: bash calcGC_umd31.sh UMD3.1_test.fa fasta_formatter  

############################################################################################
# Bovine Reference Genome is downloaded from "http://128.206.12.216/drupal/sites/bovinegenome.org/files/data/umd3.1/UMD3.1_chromosomes.fa.gz" 
# wget "URL"
# Test file preparation:
# zgrep -n '^>' UMD3.1_chromosomes.fa.gz | head -31 > umd31_rowNumbers
# awk -F':' '{print $1, $1+100}'  umd31_rowNumbers > interval_umd31_rowNumbers
# while read lines; do STARTs=`echo $lines| awk '{print $1}'`; ENDs=`echo $lines| awk '{print $2}'`; zcat UMD3.1_chromosomes.fa.gz | awk -v var1=$STARTs -v var2=$ENDs 'NR>var2 {exit} NR>=var1 && NR<=var2' >> UMD3.1_test.fa; done <interval_umd31_rowNumbers &
# bgzip UMD3.1_test.fa
#############################################################################################

umd31=${1}

# "UMD3.1_chromosomes.fa.gz" is compressed with gzip.
gzip -d ${umd31}

# To Prepare reference FASTA file with 100bp per line using "fasta_formatter" from "FASTX Toolkit"
formattedFASTA=$(basename ${umd31} .fa.gz).Formatted.fa

# Tool
fasta_formatter=${2}
${fasta_formatter} \
	-w 100 \
	-i $(basename ${umd31} .gz) \
	-o ${formattedFASTA}

# FASTA header formating
grep '^>' ${formattedFASTA} | \
	sed -e 's:^>gnl|UMD3.1|::g' -e 's:Chromosome :Chr:g' | \
	awk '{print $1"\t"$2}' \
	> fastaHeader.txt

# Calculate GC-content for bins of 100bp
outGCfile=umd31_GC_content.bed

python calcGCcontent.py \
	-i ${formattedFASTA} \
	--fastaHeader fastaHeader.txt \
	-o ${outGCfile}

########################
# Remove file
# rm ${formattedFASTA}
# Compress Fasta file 
gzip $(basename ${umd31} .gz)

######################################## END




