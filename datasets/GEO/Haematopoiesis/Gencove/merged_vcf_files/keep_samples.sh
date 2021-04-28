#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=48:00:00
#SBATCH -p shared

#INDIV="$1" snakemake -s Snakefile --cores 24


ml vcftools

vcftools --keep Gencove_samples.txt --gzvcf merged.vcf.gz --maf 0.05 --out GBR.samples.maf005 --recode

filename=GBR.samples.maf005.recode.vcf
for chromosome in {1..22}
do
	vcftools  --vcf ${filename}  --chr chr${chromosome}  --recode --recode-INFO-all --out GBR.samples.maf005.chr${chromosome}
	bgzip GBR.samples.maf005.chr${chromosome}.recode.vcf
	tabix -f GBR.samples.maf005.chr${chromosome}.recode.vcf.gz
done
