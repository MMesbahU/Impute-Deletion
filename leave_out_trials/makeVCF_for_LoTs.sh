#!/bin/bash

#
CHROM=${1}
delID=${2}
del_START=${3}
del_END=${4}
vcfFile=${5}
OutputDir=${6}
testNam=${7}
#

#
if (( ${del_START} - 1000000 < 1 )); then
        start_pos_b4_del=1
        end_pos_b4_del=$(( ${del_START} - 500 ))
        start_pos_after_del=$(( ${del_END} + 500 ))
        end_pos_after_del=$(( ${del_END} + 1000000 ))
else
        start_pos_b4_del=$(( ${del_END} - 1000000 ))
        end_pos_b4_del=$(( ${del_START} - 500 ))
        start_pos_after_del=$(( ${del_END} + 500 ))
        end_pos_after_del=$(( ${del_END} + 1000000 ))

fi

echo -e "Extracting WGS variants upstream of the deletion for Phasing \n"
#  WGS variants before deletion:
vcftools \
        --gzvcf ${vcfFile} \
        --thin 5 \
        --mac 5 \
        --recode \
        --recode-INFO-all \
        --chr ${CHROM} \
        --from-bp ${start_pos_b4_del} \
        --to-bp ${end_pos_b4_del} \
        --out ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}

echo -e "Extracting WGS variants after the deletion for phasing \n"
#  WGS variants after deletion:
vcftools \
        --gzvcf ${vcfFile} \
        --thin 5 \
        --mac 5 \
        --recode \
        --recode-INFO-all \
        --chr ${CHROM} \
        --from-bp ${start_pos_after_del} \
        --to-bp ${end_pos_after_del} \
        --out ${OutputDir}/WGS_markers_after${CHROM}_${del_START}_${del_END}
# get deletions
vcftools \
        --gzvcf ${vcfFile} \
        --snp ${delID} \
        --recode \
        --recode-INFO-all \
        --out ${OutputDir}/DEL.${delID}
#### Combine all three files
echo -e "Combining three VCF files\n"

# Combine all three files and sort
vcf-concat \
        ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}.recode.vcf \
        ${OutputDir}/DEL.${delID}.recode.vcf \
        ${OutputDir}/WGS_markers_after${CHROM}_${del_START}_${del_END}.recode.vcf | \
        vcf-sort -c \
        > ${OutputDir}/${testNam}.${delID}.vcf
# unphase
sed -i -e 's:1|1:1/1:g' -e 's:1|0:0/1:g' -e 's:0|1:0/1:g' -e 's:0|0:0/0:g' ${OutputDir}/${testNam}.${delID}.vcf
#
gzip ${OutputDir}/${testNam}.${delID}.vcf
#
rm ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}.recode.vcf ${OutputDir}/DEL.${delID}.recode.vcf ${OutputDir}/WGS_markers_after${CHROM}_${del_START}_${delEND}.recode.vcf
