#!/bin/bash

sample="$1"
cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/
VCF=${sample}.filtered.recode.vcf.gz
 
source ~/.bash_profile
ml gcc/5.5.0

refPanel=/work-zfs/abattle4/lab_data/imputation_reference_panel/10.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz
vcf_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/GRCh37

vcfFn=${vcf_dir}/${sample}.filtered.recode.GRCh37.vcf
minimac4 --refHaps ${refPanel} \
         --haps ${vcfFn} \
         --prefix testRun
