#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=48:00:00
#SBATCH -p shared

#INDIV="$1" snakemake -s Snakefile --cores 24


ml bcftools
ml vcftools

vcftools --exclude multiallelic_snps.txt --gzvcf batch4.maf005.recode.vcf.gz --out batch4.unique.samples --recode
