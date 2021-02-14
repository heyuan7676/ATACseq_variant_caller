#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH -p skylake

chromosome="$1"
python construct_VCF.py ${chromosome}

ml htslib
VCF_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/minDP3
VCFfn=${VCF_dir}/dosage_by_sample_matrix_chr${chromosome}.recode.forFastQTL.vcf
bgzip -c  ${VCFfn} > ${VCFfn}.gz
tabix ${VCFfn}.gz

VCF_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/minDP3
VCFfn=${VCF_dir}/dosage_by_sample_matrix_chr${chromosome}.recode.forFastQTL.vcf
bgzip -c  ${VCFfn} > ${VCFfn}.gz
tabix ${VCFfn}.gz

VCF_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Integration/minDP3
VCFfn=${VCF_dir}/dosage_by_sample_matrix_chr${chromosome}.recode.forFastQTL.vcf
bgzip -c  ${VCFfn} > ${VCFfn}.gz
tabix ${VCFfn}.gz
