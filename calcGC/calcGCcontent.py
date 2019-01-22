#!/usr/bin/env python

# python3.7
import os
import subprocess
import argparse
import pandas as pd

## Argument parser
parser = argparse.ArgumentParser(description="Process DNA sequences for calculating GC%.")                                               
parser.add_argument("-i", "--inputFASTA", type=str, required=True, help="Input reference sequence file in FASTA format.")
parser.add_argument("-fh", "--fastaHeader", type=str, required=True, help="Fasta file header")
parser.add_argument("-o", "--output", type=str, required=True, help="Output GC% to a file")
args = parser.parse_args()


# GC_percent function
# output: GC% and Number_of_N_bases
def GC_percent(DNA):
	N_Bases = DNA.count('n') + DNA.count('N')
	GC = float(DNA.count('c') + DNA.count('C') + DNA.count('g') + DNA.count('G')) * 100.0/(len(DNA)- N_Bases)
	return(round(GC,0), N_Bases)
##############

#### Chromosome dictionary
UMD31_header = pd.read_table(args.fastaHeader, sep='\t', header=None)

ChromDict = dict( zip(UMD31_header[0], UMD31_header[1]) )
####

# BED file Header
Proc1 = subprocess.Popen('echo -e "#Chrom\tSTART\tEND\tGC%\ttot_N_bases" > ' + args.output, shell=True)
os.waitpid(Proc1.pid, 0)

# 
f = open(args.inputFASTA, 'r+')

# Append outputs to an existing file
o = open(args.output,'a+')

lineNum=0	
	
for line in f:
	if line.startswith('>'):
		lineNum=0
		a=line.strip().split(' ')[0]
		b=a.strip().split('|')[2]
		chrom=ChromDict[b]
		continue
	else:
		gc_n=GC_percent(line)
		lineNum += 1
		o.write("\t".join([str(chrom), str((lineNum*100)-99), str(lineNum*100), str(gc_n[0]), str(gc_n[1]) ]) + "\n")
f.close()
o.close()

# END 




