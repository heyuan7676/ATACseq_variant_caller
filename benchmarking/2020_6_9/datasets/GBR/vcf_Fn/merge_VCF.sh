#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=24:00:00
#SBATCH -p shared

ml bcftools
bcftools concat ALL.chr*vcf -Oz -o ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.GBRsampleRemoved.maf005.recode.vcf.gz
