
# Set the folder to current directory 
workDIR <- system("pwd", intern=T)

setwd(workDIR)

# Rscript $outDir/../calcDPstats.R $TMPDIR/Chr${chr}.AD.FORMAT $outDir Chr${chr} 

library(data.table)

# Reading inputs: 1. ReadDepth_file 2.Output_file 
args <- commandArgs(TRUE)

#

d <- fread( cmd=paste0("zcat ", args[1]), header=T, stringsAsFactors=F, sep='\t')  

# remove columns:1-3 [Chromosome, Pos0, Pos1] 
dat <- d[,-c(1:3)]

# 
Results = t( rbind( round(apply(dat, 2, mean),2), round(apply(dat, 2, var), 2) )) 

# Columns : "SampleID", "mu_RD", "VAR_RD"

write.table(Results, file=args[2], row.names = TRUE, col.names = FALSE, quote = FALSE, sep = '\t')

