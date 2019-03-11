
## This Script will run Deletion genotyping and produce two output files
## 1st: Raw output with 13-columns; headers as: "W_0", "W_1", "W_2", "MU_0", "MU_1","MU_2","Var_0","Var_1","Var_2","ll_cur","ll_new", "No.Iter", "IndividualID"  
## 2nd: VCF file with Phred-scaled deletion genotype likelihood as "phred-likelihood of CN2,CN1,CN0" 

########################## Taking input from command line ###############################
#
commandLineArgs=commandArgs(TRUE) 
if (length(commandLineArgs) != 11){
	cat("Please provide 11 arguments after the R_script as follows:\n")
	cat("Rscript readDepth_genotyping.R <EM_script> <Read_Depth_input> <IndividualID_same_order_as_in_Read_Depth_input> <Initial_values_Mus_Vars> <Output_file> <Output_vcf> <Chromosme> <Deletion_Start> <Deletion_ID> <Reference_Allele> <Alternative_Allele>\n")
	cat(gsub(pattern = ", ", replacement = "\t", x=toString(c("Your provided inputs are: ", commandLineArgs,"\n"))))
	q(save = "no", status = 1)
}

# Under flow
min_num <- 1e-323
# Minumum variance 
min_VAR <- 1e-20

# 11 arguments supplied in "deletion_genotyping_pipeline.sh" as:
# Rscript ${ScriptDir}/readDepth_genotyping.R 11_arguments
# (1) EM script
# (2) File with Read Depth within the deletion locus for GMM input. (VCFtools's output; With header; no_of_col = 2+no_of_individual)
# (3) Individual ID as in input Read Depth file
# (4) File with Expected_RD and Variance as initial value to EM. 
# (5) Output CN class weights in a tab-separeted text file. 
# (6) Output VCF file with Phred-scaled deletion genotype likelihoods 
# (7) Chromosme 
# (8) Deletion Start Position
# (9) Deletion ID 
# (10) Reference Allele
# (11) Alternative Allele
 
##########################################################################################

# Load EM_for_GMM script  
source(commandLineArgs[1])

# Load "data.table" library
library(data.table)

############## Reading Input files ################### 

# Read-depth input:- Chrom:POS:IndividualID1...n
# Read-depth data from 3rd column ....
# we remove the header line from VCFtools output
d <- fread(commandLineArgs[2], skip = 1, stringsAsFactors = F, header = F, sep = '\t')

# Individual IDs in ReadDepth file
invIDs <- read.table(commandLineArgs[3], stringsAsFactors = F, header = F)

# Initial values for Mean and Variance: column1=IndividualID, column2=Mean_RD, column3=Variance_of_RD
expRD <- read.table(commandLineArgs[4], stringsAsFactors = F, header = F, sep='\t') 

# Expected Read Depth and Variance
expectedRD <- as.data.frame(matrix(NA, nrow = nrow(expRD), ncol = 3 ))

names(expectedRD) <- c("IndividualID","xRD","xVar")

expectedRD$IndividualID <- invIDs$V1

# Loop for filling xRD and xVar

for(i in 1:nrow(expectedRD)){
	for(j in 1:nrow(expRD)) {
		if(expectedRD$IndividualID[i] == expRD$V1[j]){ 
			expectedRD$xRD[i] <- expRD$V2[j]
			expectedRD$xVar[i] <- expRD$V3[j]
			}
	} 
}

############## Read-depth Genotyping #################

Result <- vector("list", nrow(expectedRD))

for(i in 1:nrow(expectedRD)){
  # 
  cat(paste0("EM iteration for ", expectedRD$IndividualID[i],":\n"))
  cat("No.Iter\tlog_likelihood\n")

  # Read-depth data starts from 3rd column
  dat <- d[[i+2]] 
  
  # Equal initial weight for each Copy-number (CN) [EM will update iteratively using the data]
  lambda <- rep(1/3, 3)
  
  # Expected RD for CN0:"0", CN1:"Half_of_GC_adjusted_genome-wide_average", CN2:"GC_adjusted_Genome-wide_average"
  # We will keep the expected Read-depth fixed for each CN class
  MU_expected <- c(0, expectedRD$xRD[i]/2, expectedRD$xRD[i])
  
  # Initial Variance: small positive value for CN0 [EM will update iteratively using the data]
  VAR_init <- c(0.2, expectedRD$xVar[i]/2, expectedRD$xVar[i])
  
  # EM step
  Result[[i]]  <- EM_geno(Data = dat, W_init = lambda, MU_init = MU_expected, VAR_init = VAR_init)
  cat("\n")
 
}

# Extracting Results from the list 

res <- as.data.frame(matrix(unlist(Result), nrow = nrow(expectedRD), ncol = 12, byrow = T))

names(res) <- c("W_0", "W_1", "W_2", "MU_0", "MU_1", "MU_2", "Var_0", "Var_1", "Var_2", "ll_cur", "ll_new", "No.Iter")

# Addition of corresponding individual IDs to the Result file
res$IndividualID <- expectedRD$IndividualID

# 1st Output file

write.table(res, commandLineArgs[5], col.names = T, row.names = F, quote = F, sep = '\t')

########################### Phred-scaled Genotype likelihood ############################
### 
### PL="Phred-scaled likelihoods for genotypes"
### phred score from p-value, phred = -10 * log10(p)
### probability, p = 10^-(phred/10)
### 
phred <- function(p){
  round( -10 * log10(p) ) 
}
# Maximum phred scaled likelihood 
max_PL <- phred(min_num)

#########################################################################################

no_column <- length(res$IndividualID) + 9

PLs <- as.data.frame(matrix(NA, ncol = no_column, nrow = 1))

PLs[1,c(1:9)] <- c(commandLineArgs[7], commandLineArgs[8], commandLineArgs[9], commandLineArgs[10], commandLineArgs[11], rep(".",3), "PL" )

for(i in 1:nrow(res)){
	if(res$Var_0[i] == min_VAR || res$Var_1[i] == min_VAR || res$Var_2[i] == min_VAR ){
		PLs[1, 9+i] <- "." 
	} else if(res$W_2[i] != 0 && res$W_1[i] != 0 && res$W_0[i] != 0){
		PLs[1, 9+i] <- gsub(pattern = " ", replacement = "", x = toString(phred(as.numeric(res[i, c("W_2", "W_1", "W_0")]))))
	} else if(res$W_2[i] == 0 && res$W_1[i] != 0 && res$W_0[i] != 0){ 
		PLs[1, 9+i] <-  gsub(pattern = " ", replacement = "", x = toString(c(max_PL, phred(as.numeric(res[i, c("W_1", "W_0")])) ))) 
	} else if(res$W_2[i] != 0 && res$W_1[i] == 0 && res$W_0[i] != 0){ 
		PLs[1, 9+i] <-  gsub(pattern = " ", replacement = "", x = toString(c(phred(res$W_2[i]), max_PL, phred(res$W_0[i]) ) ))
	} else if(res$W_2[i] != 0 && res$W_1[i] != 0 && res$W_0[i] == 0){
		PLs[1, 9+i] <-  gsub(pattern = " ", replacement = "", x = toString(c(phred(res$W_2[i]), phred(res$W_1[i]), max_PL ) ))
	}	else PLs[1, 9+i] <- "."
}

## Populating the VCF file 

# Adding meta-information lines 
cat("##fileformat=VCFv4.2", file = commandLineArgs[6], append = TRUE, fill = TRUE)
cat(gsub(pattern = ", ", replacement = "", x = toString( c("##fileDate=",format(Sys.Date(), "%Y%m%d")))), file = commandLineArgs[6], append = TRUE, fill = TRUE)
cat("##ALT=<ID=DEL,Description=\"Deletion\">", file = commandLineArgs[6], append = TRUE, fill = TRUE)
cat("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods (without normalization; rounded to the closest integer)\">", file = commandLineArgs[6], append = TRUE, fill = TRUE)

# Adding VCF header 
cat(gsub(pattern = ", ", replacement = "\t", x=toString(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", res$IndividualID))), file = commandLineArgs[6], append = TRUE, fill = TRUE)

# Adding Phred-scaled Deletion Genotype to VCF file
cat(gsub(pattern = ", ", replacement = "\t", x = toString( PLs[1,]) ), file = commandLineArgs[6], append = TRUE, fill = TRUE)

###################################### END ###################################################################################




