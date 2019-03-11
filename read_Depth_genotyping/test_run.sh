#!/bin/bash

# Output directory
currentDir=`pwd`
mkdir -p ${currentDir}/Results
# Tools needed
VCFTOOLS=`which vcftools`
VCF_Shuffle=`which vcf-shuffle-cols`
VCF_Concat=`which vcf-concat`
VCF_Sort=`which vcf-sort`
RSCRIPT=`which Rscript`
JAVA=`which java`
BEAGLE=b4.r1274.jar
PERL=`which perl`
outDir=${currentDir}/Results
scriptDir=ScriptDir
VCFfile=test.Chr21.vcf.gz
expected_read_depth=Expected_Read_Depth.txt
# Run the script 
bash ${scriptDir}/deletion_genotyping_pipeline.sh ${scriptDir} ${outDir} ${VCFfile} Chr21 21184869 21188202 C "<DEL>" ${expected_read_depth} ${VCFTOOLS} ${VCF_Shuffle} ${VCF_Concat} ${VCF_Sort} ${RSCRIPT} ${JAVA} ${BEAGLE} ${PERL} 1>>test.log 2>>test.err

