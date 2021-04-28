#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=48:00:00
#SBATCH -p shared

#INDIV="$1" snakemake -s Snakefile --cores 24


ml bcftools
ml vcftools

#bcftools query -l batch2.maf005.biallelic.recode.vcf.gz | sort > batch2_samples
#bcftools query -l batch3.maf005.biallelic.recode.vcf.gz | sort > batch3_samples

#grep -F -f batch2_samples batch3_samples > duplicated_samples.txt
vcftools --remove  duplicated_samples.txt --gzvcf batch4.maf005.recode.vcf.gz --out batch4.unique.samples --recode
