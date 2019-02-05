#!/bin/bash

echo -e "Started at $(date) \n\n"

# To run in Slurm for all chromosome 
# for chr in {1..29}; do sbatch run_extractReadDepth.sh outDir Chr${chr}.vcffile GC_annotation samples_2_keep ${chr}; done
outDir=${1}
VCFfile=${2}
GC_annotation=${3}
samples_2_keep=${4}
chr=${5}

mkdir -p ${outDir}/tmpRD

echo -e "Provided files: ${VCFfile} : ${GC_annotation} : ${samples_2_keep} \n"

echo -e "Extract read depth from chromosome ${chr} \n"

# Extract read depth from each chromosome
echo -e "output in ${outDir}/tmpRD directory \n"

bash extractReadDepth.sh ${VCFfile} ${GC_annotation} ${samples_2_keep} ${chr} ${outDir}/tmpRD

echo -e "Ended at $(date) \n\n"


