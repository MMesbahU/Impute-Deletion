#!/bin/bash

# Summary stats for 29 autosomes for all samples
# bash run_calcRDstats.sh < 1. GC annotation bed file from step 1.2 > < 2. output directory from step 1.3 > < 3. bedfile_header file from step 1.3 >
# expected inputs: 1. GC_annotation file 2. directory of RD output file 3. bedfile_header
GC_annotation=${1}
dirGCfiles=${2}
# from  ${dirGCfiles}/Chr${chr}.bedfile_header.txt 
bedheader=${3}
# Calculate mean and variance per GC%

while read lines

do 
	awk -v OFS='\t' 'NR==1 {print $0}' ${bedheader} > ${dirGCfiles}/ReadDepth.Chr1_29.GC${lines}.txt
	
	numCols=`awk 'NR==1 {print NF}' ${dirGCfiles}/ReadDepth.Chr1_29.GC${lines}.dat`

	# Removing redundant column 1-3 	
	cut -f4-${numCols} ${dirGCfiles}/ReadDepth.Chr1_29.GC${lines}.dat >> ${dirGCfiles}/ReadDepth.Chr1_29.GC${lines}.txt

	gzip ${dirGCfiles}/ReadDepth.Chr1_29.GC${lines}.txt

	rm ${dirGCfiles}/ReadDepth.Chr1_29.GC${lines}.dat

# output file: column1=IndividualID(as in VCF) "\t" Column2=Average_ReadDepth_for_GC% "\t" Column3=Avg_Variance_for_GC% #
	Rscript calcRDstats.R ${dirGCfiles}/ReadDepth.Chr1_29.GC${lines}.txt.gz ${dirGCfiles}/Expected_RD_GC${lines}.txt 

done < <(awk '{print $NF}' ${GC_annotation} | sort -V | uniq )


