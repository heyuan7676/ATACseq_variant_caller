#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=48:00:00
#SBATCH -p lrgmem

#INDIV="$1" snakemake -s Snakefile --cores 24


ml bcftools

filename=merged.bcf
bcftools +fill-tags  ${filename} -- -t AF | grep -v "#" | awk '{print $1,$2,$8}' > ${filename/.vcf.gz/}.maf.txt
