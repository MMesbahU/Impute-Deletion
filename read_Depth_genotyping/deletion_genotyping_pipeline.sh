#!/bin/bash

## Clock timing
echo -e "\n\n"
echo -e "Job started at: $(date) \n\n"
Job_START=$(date +%s)

################################ Inputs expected #############################################################################
# (1) Script directory for R, Beagle, etc, e.g. ~/deletion_genotyping/Scripts												 #
# (2) Output directory, e.g. ~/deletion_genotyping/Outputs                                                                   #
# (3) VCF file with location, e.g. ~/deletion_genotyping/VCFdir/vcffile.vcf.gz                                               #
# (4) Chromosome, e.g Chr1 or 1 (as in VCF file)                                                                             #
# (5) Deletion START position, e.g. 12345670                                                                                 #
# (6) Deletion END position, e.g. 12355670                                                                                   #
# (7) Reference Allele, e.g. A                                                                                               #
# (8) Alternative allele, e.g. <DEL>                                                                                         #
# (9) File containing Expected read-depth and variance for diploid genome of each individual (for certain GC% content), e.g. #
#	~ No header is needed                                                                                                    #
#	~ One Row for each individual                                                                                            #
#	~ Three columns with tab separation ("\t"):                                                                              #
#		* column1=IndividualID (as in VCF)                                                                                   #
#		* Column2=Average ReadDepth for individual for certain GC-content (see GC_calculation_from_VCF_pipeline.sh)          #
#		* Column3=Avg. Variance for individual for certain GC-content		                                                 #
# Followed by Eight softwares with full path.                                                                                #
# To run the script                                                                                                          #
# bash deletion_genotyping_pipeline.sh \                                                                                     #
#	~/deletion_genotyping/Scripts \                                                                                          #
#	~/deletion_genotyping/Outputs \                                                                                          #
#	~/deletion_genotyping/VCFdir/vcffile.vcf.gz \                                                                            #
#	Chr1 \                                                                                                                   #
#	12345670 \                                                                                                               #
#	12355670 \                                                                                                               #
#	A \                                                                                                                      # 
#	<DEL> \                                                                                                                  # 
#	~/FilePath/ExpectedReadDepth.txt \                                                                                       #
#	(10) ... (17) \                                                                                                          #
#	1>>logFile 2>>errFile                                                                                                    # 
############################################################################################################################## 

#################### Checking inputs 
if [[ ${#*} != 17 ]]; then
	echo -e "Please provide 17 arguments after bash_script as described below:"
	echo -e "\t(1) Script directory for R, Beagle, etc, e.g. ~/deletion_genotyping/Scripts"
	echo -e "\t(2) Output directory, e.g. ~/deletion_genotyping/Outputs"
	echo -e "\t(3) VCF file with location, e.g. ~/deletion_genotyping/VCFdir/vcffile.vcf.gz"
	echo -e "\t(4) Chromosome, e.g Chr1 or 1 (as in VCF file)"
	echo -e "\t(5) Deletion START position, e.g. 12345670"
	echo -e "\t(6) Deletion END position, e.g. 12355670"
	echo -e "\t(7) Reference Allele, e.g. A"
	echo -e "\t(8) Alternative allele, e.g. <DEL>"
	echo -e "\t(9) File containing Expected read-depth and variance for diploid genome of each individual (for certain GC% content)" 
	echo -e "\t In addition, please provide full path of follwoing softwares (or edit the script accordingly):"
	echo -e "\t\t (10) vcftools (11) vcf-shuffle-cols (12) vcf-concat (13) vcf-sort (14) Rscript (15) java (16) beagle (17) perl"
	exit 1
fi
###########################
#### Bash Inputs
ScriptDir=${1}
OutputDir=${2}
vcfFile=${3}
CHROM=${4}
del_START=${5}
del_END=${6}
RefAllele=${7}
AltAllele=${8}
ExpRDfile=${9}
VCFTOOLS=${10}
VCF_Shuffle=${11}
VCF_Concat=${12}
VCF_Sort=${13}
RSCRIPT=${14}
JAVA=${15}
BEAGLE=${16}
PERL=${17}
delID=${CHROM}:${del_START}_${del_END}

#
echo -e "Your inputs are:\n \t(1) ${ScriptDir}\n\t(2) ${OutputDir}\n\t(3) ${vcfFile}\n\t(4-6) ${delID}\n\t (7-8) ${RefAllele}|${AltAllele}\n\t (9) ${ExpRDfile}\n\t (10-17) and eight Softwares with full path!! \n\n"

# Individual ID for analysis
awk '{print $1}' ${ExpRDfile} > ${OutputDir}/Individual_IDs.list
InvIDs=${OutputDir}/Individual_IDs.list

##############################################################################################
###############                  (._.)  Extract read-depth from VCF file  (._.)       ########
###############         				for the Deletion region                       ########
############### 																	  ########
###############				 					      								  ########
##############################################################################################
echo -e "Extracting read-depth within deletion ${delID} from the VCF file ${vcfFile} \n"
ReadDepth4GMM_START=$(date +%s)
### VCFTOOLS to extract Allelic-Depth 
# Read depth info can be extracted from DP or AD from FORMAT column. 
# We used AD-tag from the 1000 Bull Genomes Project VCF file
# output file will have file extension ".AD.FORMAT"
${VCFTOOLS} \
	--gzvcf ${vcfFile} \
	--chr ${CHROM} \
	--from-bp ${del_START} \
	--to-bp ${del_END} \
	--minQ 30 \
	--min-alleles 2 \
	--max-alleles 2 \
	--remove-indels \
	--thin 10 \
	--keep ${InvIDs} \
	--extract-FORMAT-info AD \
	--out ${OutputDir}/ReadDepth_${CHROM}_${del_START}_${del_END}

ReadDepth4GMM_END=$(date +%s)
# This step is for adding ADF and ADR. (Not needed if DP information is used instead)
# Adding forward and reverse reads (ADF="Allelic depths on the forward strand" and ADR="Allelic depths on the reverse strand") [Perl code credit: Prof. Bernt Guldbrandtsen)
cat ${OutputDir}/ReadDepth_${CHROM}_${del_START}_${del_END}.AD.FORMAT | \
	${PERL} -ane 'if ($.==1) {print join("\t",@F),"\n"; next} print shift @F,"\t",shift @F;while(@F){ $v=shift @F;($v1,$v2)=split(",",$v); print "\t",$v1+$v2} print "\n";' \
	> ${OutputDir}/ReadDepth_${CHROM}_${del_START}_${del_END}.txt
##############################################################################################

####################################################################################################
###############          (._.)  Estimating Copy-Number Likelihood  (._.)               		########
###############          using Expectation Maximization (EM) algorithm                		########
###############        constrained Gaussian Mixture model (with fixed mean)    			  	########
###############        For a Deletion locus: three copy-number classes CN0, CN1, and CN2    ########
####################################################################################################
echo -e "Estimating Copy-Number Likelihood using constrained Gaussian Mixture model\n\twith three copy-number (CN) classes: CN0, CN1, CN2 \n\n"
# Read depth file for genotyping		
readDepth=${OutputDir}/ReadDepth_${CHROM}_${del_START}_${del_END}.txt
# Individual IDs as in input Read-depth 
awk 'NR==1' ${readDepth} | sed 's:\t:\n:g' | sed '1,2d' > ${OutputDir}/inv_list_fromReadDepth_${CHROM}_${del_START}_${del_END}.list
inv_list_fromRD=${OutputDir}/inv_list_fromReadDepth_${CHROM}_${del_START}_${del_END}.list
# output files from R 
ouput_Deletion_Weights=${OutputDir}/DelWeight_${CHROM}_${del_START}_${del_END}.txt
vcf_Phred_likelihood=${OutputDir}/PL_${CHROM}_${del_START}_${del_END}.vcf
##Running R script for GMM
GMM_script=${ScriptDir}/EM_GMM.R

# Timming
ReadDepth_genotyping_START=$(date +%s)
# Running R for genotyping
${RSCRIPT} ${ScriptDir}/readDepth_genotyping.R \
	${GMM_script} \
	${readDepth} \
	${inv_list_fromRD} \
	${ExpRDfile} \
	${ouput_Deletion_Weights} \
	${vcf_Phred_likelihood} \
	${CHROM} \
	${del_START} \
	${delID} \
	${RefAllele} \
	${AltAllele}
#
ReadDepth_genotyping_END=$(date +%s)
#################################################### GMM END ###########################################

######################################################################################################################
###############    (._.)  Leveraging flanking SNP haplotypes to improve Deletion Genotypes  (._.)             ########
###############                 Phred-scaled deletion likelihood (PL) to integer Genotype                     ########
###############                           using BEAGLE software (V4)                                          ########
###############                                                                                               ########
######################################################################################################################

###### Filtering WGS Variants
# Keep SNPs/indels when: 
#	-maximum missing 10% or less  
#	--minQ 30:  Include only sites with Phred Quality value above 30
#	--thin 5: Thining of 5bp; i.e. no two variants are within 5 bp from one another
# 	--mac 5: minimum allele count is 5 (and more)
#	--max-missing 0.1: for a marker, a maximum of 10% (or less) individuals missing genotypes
#	--min & max-alleles 2: Only bi-allelic variants
########
# Extract biallelic WGS variants within 1Mb of putative Start/END of the deletion,
#	 	after excluding markers within 500bp of the deletion breakpoint (to account for breakpoint uncertainty)
start_pos_b4_del=$(( ${del_START} - 1000000 ))
end_pos_b4_del=$(( ${del_START} - 500 ))
start_pos_after_del=$(( ${del_END} + 500 ))
end_pos_after_del=$(( ${del_END} + 1000000 ))
#
SNPextract4Beagle_START=$(date +%s)
echo -e "Extracting WGS variants upstream of the deletion for Phasing \n"
#  WGS variants before deletion: 
${VCFTOOLS} \
	--gzvcf ${vcfFile} \
	--minQ 30 \
	--thin 5 \
	--mac 5 \
	--max-missing 0.1 \
	--min-alleles 2 \
	--max-alleles 2 \
	--keep ${inv_list_fromRD} \
	--recode \
	--recode-INFO-all \
	--chr ${CHROM} \
	--from-bp ${start_pos_b4_del} \
	--to-bp ${end_pos_b4_del} \
	--out ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}

##### HWE filtering: ######################
# vcftools --vcf ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}.recode.vcf --hardy --out ${OutputDir}/HWE_WGS_markers_before${CHROM}_${del_START}_${del_END}
# awk 'NR>1 && $NF<1e-5 {print $1"\t"$2}' ${OutputDir}/HWE_WGS_markers_before${CHROM}_${del_START}_${del_END}.hwe > ${OutputDir}/SNPs_2_remove_b4_DEL.txt
# vcftools --vcf ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}.recode.vcf --exclude-positions ${OutputDir}/SNPs_2_remove_b4_DEL.txt --recode --recode-INFO-all --out ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}_HWE
##########################################
echo -e "Extracting WGS variants after the deletion for phasing \n"
#  WGS variants after deletion: 
${VCFTOOLS} \
	--gzvcf ${vcfFile} \
	--minQ 30 \
	--thin 5 \
	--mac 5 \
	--max-missing 0.1 \
	--min-alleles 2 \
	--max-alleles 2 \
	--keep ${inv_list_fromRD} \
	--recode \
	--recode-INFO-all \
	--chr ${CHROM} \
	--from-bp ${start_pos_after_del} \
	--to-bp ${end_pos_after_del} \
	--out ${OutputDir}/WGS_markers_after${CHROM}_${del_START}_${del_END}
#
SNPextract4Beagle_END=$(date +%s)

#### Combine all three files 
echo -e "Combining three VCF files\n"
# Ordering VCF Columns
${VCF_Shuffle} \
	--template ${vcf_Phred_likelihood} \
	${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}.recode.vcf \
	> ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}.vcf
#
${VCF_Shuffle} \
	--template ${vcf_Phred_likelihood} \
	${OutputDir}/WGS_markers_after${CHROM}_${del_START}_${del_END}.recode.vcf \
	> ${OutputDir}/WGS_markers_after${CHROM}_${del_START}_${del_END}.vcf
# Combine all three files and sort 
${VCF_Concat} ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}.vcf ${vcf_Phred_likelihood} ${OutputDir}/WGS_markers_after${CHROM}_${del_START}_${del_END}.vcf | \
	${VCF_Sort} -c \
	> ${OutputDir}/WGS_markers_for_beagle${CHROM}_${del_START}_${del_END}.vcf


########## Refine Deletion Genotype using Beagle ###########
echo -e "Refining Deletion Genotype using Beagle\n"
beagle_input=${OutputDir}/WGS_markers_for_beagle${CHROM}_${del_START}_${del_END}.vcf

beagle_output=${OutputDir}/Phased_WGS_markers_for_beagle${CHROM}_${del_START}_${del_END}

Beagle_START=$(date +%s)

${JAVA} -Xmx1024m -jar ${BEAGLE} \
	gl=${beagle_input} \
	out=${beagle_output} \
	burnin-its=10 phase-its=15 window=12000 overlap=2000 nthreads=1

Beagle_END=$(date +%s)
################# Phasing ENDs  ###################################

########################################################################################################### 
# 			Testing Del genotype accuracy without data augmentation using EM_geno step:
#
# cat ${OutputDir}/WGS_markers_before${CHROM}_${del_START}_${del_END}.vcf > ${OutputDir}/Without_EM_step.WGS_markers_for_beagle${CHROM}_${del_START}_${del_END}.vcf
# indvNum=`wc -l ${inv_list_fromRD} | awk '{print $1}'`
# paste <(grep -w ${delID} ${vcf_Phred_likelihood} | cut -f1-9) \
#       <(yes . | head -${indvNum} | tr '\n' '\t') \
#       >> ${OutputDir}/Without_EM_step.WGS_markers_for_beagle${CHROM}_${del_START}_${del_END}.vcf
# grep -v '^#' ${OutputDir}/WGS_markers_after${CHROM}_${del_START}_${del_END}.vcf >> ${OutputDir}/Without_EM_step.WGS_markers_for_beagle${CHROM}_${del_START}_${del_END}.vcf
# Beagle 
# ${JAVA} -Xmx1024m -jar ${BEAGLE} gl=${OutputDir}/Without_EM_step.WGS_markers_for_beagle${CHROM}_${del_START}_${del_END}.vcf out=${OutputDir}/Phased.Without_EM_step.WGS_markers_for_beagle${CHROM}_${del_START}_${del_END} burnin-its=10 phase-its=15 window=12000 overlap=2000 nthreads=1
###########################################################################################################

######## Removing intermediate files
rm -rf ${OutputDir}/ReadDepth_*
rm -rf ${OutputDir}/*.list
rm -rf ${OutputDir}/*.vcf

################################## Clock time #############################################################
Job_END=$(date +%s)
echo "Computation time details:"
echo $((ReadDepth4GMM_END - ReadDepth4GMM_START)) | awk '{print "Total time for Read Depth extraction for GMM: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
echo $((ReadDepth_genotyping_END - ReadDepth_genotyping_START)) | awk '{print "Total time for Read Depth genotyping using GMM: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
echo $((SNPextract4Beagle_END - SNPextract4Beagle_START)) | awk '{print "Total time for SNP extraction for imputation (using VCFtools): " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
echo $((Beagle_END - Beagle_START)) | awk '{print "Total time for Beagle: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
echo $(( Job_END - Job_START)) | awk '{print "Total run time for the deletion genotyping pipeline: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
echo -e "\n\n"
echo -e "Job ended at: $(date) \n\n"
###########################################################################################################


