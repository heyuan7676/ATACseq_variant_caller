#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=48:00:00
#SBATCH -p skylake

#INDIV="$1" snakemake -s Snakefile --cores 24


ml htslib
ml vcftools

chromosome="$1"
#filename=batch4.unique.samples.recode.vcf
filename=merged.bcf.gz
vcftools  --gzvcf ${filename}  --chr chr${chromosome}  --recode --recode-INFO-all --out merged.chr${chromosome}
bgzip merged.chr${chromosome}.recode.vcf
tabix -f merged.chr${chromosome}.recode.vcf.gz

