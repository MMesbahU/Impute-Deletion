#!/bin/bash

# to run 
# 1. VCF file 
# 2. GC_annotation.bed
# 3. Sample_IDs list to keep
# 4. Chr
# 5. outputDir
# bash extractReadDepth.sh VCF_file GC_annotation Sample_IDs_2_keep Chr outputDir

VCFfile=${1}
GC_annotation=${2}
samples_2_keep=${3}
chr=${4}
outputDir=${5}

# Extract Read depth from VCF file
# Filter: Bi-allelic SNPs only; SNP quality >=30; no SNPs within 10 bp of each other
vcftools \
	--gzvcf ${VCFfile} \
	--minQ 30 \
	--min-alleles 2 \
	--max-alleles 2 \
	--remove-indels \
	--thin 10 \
	--keep ${samples_2_keep} \
	--extract-FORMAT-info AD \
	--out ${outputDir}/Chr${chr}.raw
# 
cat ${outputDir}/Chr${chr}.raw.AD.FORMAT | \
	perl -ane 'if ($.==1) {print join("\t",@F),"\n"; next} print shift @F,"\t",shift @F;while(@F){ $v=shift @F;($v1,$v2)=split(",",$v); print "\t",$v1+$v2} print "\n";' > ${outputDir}/Chr${chr}.raw.AD.txt
#
totalCols=`awk 'NR==2 {print NF}' ${outputDir}/Chr${chr}.raw.AD.FORMAT`
# Header for BED file
paste -d '\t' <(echo -e "#CHROM\tPOS_0\tPOS_1") <( head -1 ${outputDir}/Chr${chr}.raw.AD.FORMAT | cut -f3-${totalCols}) > ${outputDir}/Chr${chr}.AD.bed

# for Sanity check 
cp ${outputDir}/Chr${chr}.AD.bed ${outputDir}/Chr${chr}.bedfile_header.txt
#
paste -d '\t' <(awk -v OFS='\t' 'NR>1 {print $1, $2-1, $2}' ${outputDir}/Chr${chr}.raw.AD.txt) <(awk -v OFS='\t' 'NR>1 {print $0}' ${outputDir}/Chr${chr}.raw.AD.txt | cut -f3-${totalCols} ) >> ${outputDir}/Chr${chr}.AD.bed 

gzip ${outputDir}/Chr${chr}.AD.bed

# remove 
rm ${outputDir}/Chr${chr}.raw.AD.FORMAT 
rm ${outputDir}/Chr${chr}.raw.AD.txt

# Bedtools intersect 
while read lines
do 
	intersectBed \
		-a <(awk -v OFS='\t' -v CHR=Chr${chr} -v GC=${lines} '$1==CHR && $4==GC {print $1, $2, $3}' ${GC_annotation} ) \
		-b ${outputDir}/Chr${chr}.AD.bed.gz \
		-wb >> ${outputDir}/ReadDepth.Chr1_29.GC${lines}.dat

done < <(awk '{print $NF}' ${GC_annotation} | sort -V | uniq )

# END


