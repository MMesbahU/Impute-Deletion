#!/bin/bash

# To Run: bash calcGC_umd31.sh FASTAfile.fa.gz fasta_formatter  

# Bovine Reference Genome is downloaded from "http://128.206.12.216/drupal/sites/bovinegenome.org/files/data/umd3.1/UMD3.1_chromosomes.fa.gz" 
# wget "URL"

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

# END




