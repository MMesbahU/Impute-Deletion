#!/bin/bash

# bash filter_GC_interval.sh umd31_GC_content.bed nestedRepeats.txt.gz CNV_SV_file_from_DGVA_22Jan2019.list

GCfile=${1}
repeatRegions=${2}
CNVs_study_list=${3}
##### Filtering the BED interval: 
## (1) Keeping Bovine Autosomes: Chr1-Chr29
## (2) Exclude last line from each Chromosome (it is usually < 100bp)
## (3) Exclude intervals that contain "N" bases. These are assembly gaps.

for chr in {1..29}
do 
	grep -w Chr${chr} ${GCfile} | \
	sed '$ d' | \
	awk -v OFS='\t' '$NF==0 {print $1, $2, $3, $4}' \
	>>  Chr1_29.noGAPs.umd31.bed 
done

## (4) Exclude Repeat regions 
# wget http://hgdownload.soe.ucsc.edu/goldenPath/bosTau6/database/nestedRepeats.txt.gz
# repeatRegions=nestedRepeats.txt.gz

for chr in {1..29}
do 
	zgrep -w chr${chr} ${repeatRegions} | \
	awk -v OFS='\t' '{gsub(/chr/, "Chr",$2)}{print $2, $3, $4}' | \
	sort -V | \
	uniq >>  Chr1_29_NestedRepeats.bed
done

# 
intersectBed \
	-a Chr1_29.noGAPs.umd31.bed \
	-b Chr1_29_NestedRepeats.bed \
	-v > Chr1_29.noGAPs.noRepeats.umd31.bed

## (5) Exclude known CNV and SV regions, such as from DGVA studies with UMD3.1 mapping
# the "list" includes 8 studies from DGVA. Alternatively, SVs from ENSEMBL 94 with Ref UMD3.1 
# wget ftp://ftp.ensembl.org/pub/release-94/variation/vcf/bos_taurus/bos_taurus_structural_variations.vcf.gz 
# vcftools --gzvcf bos_taurus_structural_variations.vcf.gz --get-INFO END -c | awk -v OFS='\t' 'NR>1 && $1~/^[0-9*]/ {print "Chr"$1, $2, $5,($5-$2+1)}' | awk -v OFS='\t' '$4<=1e6 {print $1, $2, $3}' | sort -V | uniq > Chr1_29.ensembl94.bos_taurus_structural_variations.bed

wget -i ${CNVs_study_list}

# Keeping CNVs or SV region <= 1MB size, in authosomes
for file in $(ls *Bos_taurus_UMD_3.1*) 
do  
	zgrep -v '^#' ${file} | \
	awk -v OFS='\t' '{print $1,$4,$5,($5-$4+1)}' | \
	sed -e 's:chr::g' -e 's:Chr::g' |  \
	awk awk -v OFS='\t' '$1~/^[0-9*]/ && $4<=1e6 {print "Chr"$1, $2, $3}' | \
	sort -V | uniq >> Chr1_29.UMD31_DGVA_CNVs_SVs.bed
done

# Excluding CNV/SV regions
intersectBed \
	-a Chr1_29.noGAPs.noRepeats.umd31.bed \
	-b Chr1_29.UMD31_DGVA_CNVs_SVs.bed \
	-v > Chr1_29.noGAPs.noRepeats.noCNVs.umd31.bed

###########################################################################
# Remove few files
# rm ${GCfile} Chr1_29.noGAPs.umd31.bed Chr1_29.noGAPs.noRepeats.umd31.bed  
# rm *Bos_taurus_UMD_3.1*.gvf.gz
##########################################################################
