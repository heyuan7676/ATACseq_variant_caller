#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=48:00:00
#SBATCH -p lrgmem

#INDIV="$1" snakemake -s Snakefile --cores 24

chr="$1"
x=chr${chr}.maf005.biallelic.recode.vcf.gz
echo $x 
zcat ${x} | grep -v "#"  | awk '{print $3, $1"_"$2}' | sed "s/chr//g" > ${x/.maf005.biallelic.recode.vcf.gz/.variants.txt}
