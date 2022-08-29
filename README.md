# **Impute large chromosomal deletions** 

## Building Reference Haplotype Panel
We used following pipelines:

### 1. Calculate GC-content in bins of 100bp for Bovine Genome Assembly UMD3.1. 

The assembly can be downloaded from [UMD3.1](http://bovinegenome.elsiklab.missouri.edu/node/61). 

1.1 Output will be similar to files in "calcGC/" directory

```
bash calcGC/calcGC_umd31.sh UMD3.1_test.fa pasth_fasta_formatter_software
```

1.2 Filter BED intervals, e.g. assembly GAPs, [repeat regions from "nestedRepeats.txt.gz"](http://hgdownload.soe.ucsc.edu/goldenPath/bosTau6/database), known cattle CNV/SV regions from [DGVA](https://www.ebi.ac.uk/dgva/data-download) or [ENSEMBL](ftp://ftp.ensembl.org/pub/release-94/variation/vcf/bos_taurus/bos_taurus_structural_variations.vcf.gz), sex chromosomes, unplaced contigs, Mitochondrial DNA etc.  

```
bash calcGC/filter_GC_interval.sh umd31_GC_content.bed nestedRepeats.txt.gz CNV_SV_file_from_DGVA_22Jan2019.list    
```
#### Calculate genome-wide avg. read depth for each animal for each GC-bin. 
1.3 Extract read depth from VCF file for each animals, for all 29 bovine autosomes. Use 'bedtools intersect' to extract overlapping GC intervals. 

```
# expected inputs: 1.output directory; 2. Chr_num_VCF_file 3.GC_annotation (bed file from step 1.2); 4. list of samples to keep from the vcf file; 5. chr number  
# for chr in {1..29}; do bash calcGC/run_extractReadDepth.sh outputDir Chr${chr}.vcf.gz GC_annotation samples_2_keep ${chr} 1>>Chr${chr}.log 2>>Chr${chr}.err; done
bash calcGC/run_extractReadDepth.sh outputDir Chr${chr}.vcf.gz GC_annotation samples_2_keep ${chr} 1>>Chr${chr}.log 2>>Chr${chr}.err
```

1.4 Calculate summary Read Depth stats for all 29 autosomes

```
bash calcGC/run_calcRDstats.sh <1.GC_annotation_bed_file_from_step_1.2> <2.output_directory_from_step_1.3> <3.bedfile_header_file_from_step_1.3> 
```
## 2. Building Haplotype Reference for imputing large chromosomal deletions

### Extend reference panel using read depth from VCF file [Read Depth genotyping]( http://pure.au.dk/portal/en/publications/genotype-call-for-chromosomal-deletions-using-readdepth-from-whole-genome-sequence-variants-in-cattle(a42d451c-ebbe-49ca-8dc0-61c166bb120c).html )




## 3. Imputation

##### References
1. Mesbah-Uddin, M., Guldbrandtsen, B., Iso-Touru, T., Vilkki, J., De Koning, D. J., Boichard, D., . . . Sahana, G. (2018). Genome-wide mapping of large deletions and their population-genetic properties in dairy cattle. DNA Res, 25(1), 49-59. (https://doi.org/10.1093/dnares/dsx037)

2. Mesbah-Uddin, M., Guldbrandtsen, B., Lund, M. S., & Sahana, G. (2018). Genotype call for chromosomal deletions using read-depth from whole genome sequence variants in cattle. Paper presented at the Proc. World Congr. Genet. Appl. Livest. Prod., Auckland, New Zealand. http://www.wcgalp.org/system/files/proceedings/2018/genotype-call-chromosomal-deletions-using-read-depth-whole-genome-sequence-variants-cattle.pdf

3. Mesbah-Uddin, M., Guldbrandtsen, B., Lund, M. S., Boichard, D., & Sahana, G. (2019). Joint imputation of whole-genome sequence variants and large chromosomal deletions in cattle. J Dairy Sci, 102(12), 11193-11206.(https://doi.org/10.3168/jds.2019-16946)

