#!/bin/bash

# first run "makeVCF_for_LoTs.sh" script; then this script 
# Run variables: bash run_leave_out_trials.sh outDir delID inputVCF maskVCF_for_LoTs.R
outDir=${1}
delID=${2}
inputVCF=${3}
maskVCF=${4}
Beagle=b4.r1274.jar
#
zgrep -vE '^#' ${inputVCF} > ${outDir}/temp.$(basename ${inputVCF} .vcf.gz).txt

gzip ${outDir}/temp.$(basename ${inputVCF} .vcf.gz).txt

input_vcf=${outDir}/temp.$(basename ${inputVCF} .vcf.gz).txt.gz

del_row_no=`zgrep -vE '^#' ${input_vcf} | grep -wn ${delID} | cut -d ':' -f1`

# 772 animals
# Run for 77 trials
echo -e "Leave-out-trials for ${delID} started at $(date) \n"

for i in {1..77}
do
    echo -e "Trial No. ${i} Started at $(date) \n"

    zgrep -E '^#' ${inputVCF} > ${outDir}/trial${i}.$(basename ${inputVCF} .vcf.gz).vcf

    outVCF=${outDir}/trial${i}.$(basename ${inputVCF} .vcf.gz).vcf

    trial_no=${i}

    output_mask=${outDir}/maskedVCF${i}.$(basename ${inputVCF} .vcf.gz).txt

    Rscript ${maskVCF} ${input_vcf} ${del_row_no} ${trial_no} ${output_mask}

    cat ${output_mask} >> ${outVCF}

    gzip ${outVCF}

    rm ${output_mask}

    ################# Phasing with Beagle  ###################################
    beagleInput=${outVCF}.gz
    beagleOutput=${outDir}/imputed.$(basename ${beagleInput} .vcf.gz)
    java -Xmx30g -jar ${Beagle} \
        gt=${beagleInput} \
        out=${beagleOutput} \
        burnin-its=10 \
        phase-its=15 \
        window=12000 \
        overlap=2000 \
        nthreads=5
    ########################################################
    rm ${beagleInput}

    echo -e "Trials NO. ${i} Finished at $(date) \n"

done

echo -e "Total 77 trials for ${delID} Finished at $(date)"

