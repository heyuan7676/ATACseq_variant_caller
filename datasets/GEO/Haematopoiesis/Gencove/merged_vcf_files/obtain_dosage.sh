#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH -p skylake

#INDIV="$1" snakemake -s Snakefile --cores 24


ml bcftools

chromosome="$1"

filename=chr${chromosome}.maf005.biallelic.recode.vcf.gz
output=chr${chromosome}.maf005.biallelic.recode.dosage.txt
bcftools query -f '%CHROM\t%POS[\t%DS]\n' ${filename} | sed 's/chr//g' | sed -r 's/\s+/_/'  > ${output}


