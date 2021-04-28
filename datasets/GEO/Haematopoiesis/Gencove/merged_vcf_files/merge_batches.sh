#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=48:00:00
#SBATCH -p skylake

#INDIV="$1" snakemake -s Snakefile --cores 24


ml htslib
ml bcftools

chromosome="$1"

filename=chr${chromosome}.maf005.vcf
bcftools merge batch4.maf005.chr${chromosome}.recode.vcf.gz merged.chr${chromosome}.recode.vcf.gz | grep -v "\./\." > ${filename}
bgzip ${filename}
tabix -f ${filename}.gz


# remove multi-allelic snps
bcftools +fill-tags  ${filename}.gz -- -t AF | grep -v "#" | awk '{print $1,$2,$3,$8}' > ${filename/.vcf/}.maf.txt
awk -F',' '{if(NF==2) print $0}'  ${filename/.vcf/}.maf.txt | awk -F' ' '{print $3}' > multiallelic_snps_chr${chromosome}.txt


#ml vcftools
#vcftools --exclude multiallelic_snps_chr${chromosome}.txt --gzvcf ${filename}.gz --out chr${chromosome}.maf005.biallelic --recode
#filename=chr${chromosome}.maf005.biallelic.recode.vcf
#bgzip ${filename}
#tabix -f ${filename}.gz


