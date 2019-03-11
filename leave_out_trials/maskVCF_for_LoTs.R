rm(list = ls())
# R version 3.5.2 (2018-12-20)
# 1.input_vcf 2. del_row_no 3.trial_no 4.output_mask
commandLineArgs <- commandArgs(TRUE)

input_vcf <- commandLineArgs[1]

del_row_no <- as.integer(commandLineArgs[2])

trial_no <- as.integer(commandLineArgs[3])

masked_VCF <- commandLineArgs[4]

# 772 animals
# indexing animals to mask the deletion genotype
# random sample 10 animals for 76 times;
# and 12 animals for last iteration (77th)
# Seed = 2019
columNumbers <- 10:781
columnIndex <- vector("list",77)
for(i in 1:length(columnIndex)){
  if( length(columNumbers)>12 ){
    set.seed(2019)
    a <- sample(columNumbers,10, replace = FALSE)
    columnIndex[[i]] <- sort(a)
    columNumbers <- columNumbers[-match(a,columNumbers)]
  } else
    columnIndex[[i]] <- columNumbers
 }
#####################################################################
# Read genotype file
# We randomly selected a group of 10 animals (out of 772 animals)
# and set all SV markers missing ="./."
# Mask genotypes
######################################################################

library(data.table)
# have to remove "cmd" data.table version  before 1.11.8
# d <- fread(cmd = paste0("zcat ", input_vcf), stringsAsFactors = F, header = F, sep='\t')
d <- fread(paste0("zcat ", input_vcf), stringsAsFactors=F, header=F, sep='\t')

d[del_row_no, columnIndex[[trial_no]]] <- "./."

write.table(d, masked_VCF, col.names = F, row.names = F, quote = F, sep='\t')
