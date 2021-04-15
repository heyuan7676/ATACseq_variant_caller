#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=24:00:00
#SBATCH -p lrgmem

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis//ATAC_seq/alignment_bowtie

ml fastqtl
ml bcftools
ml htslib
ml python/3.7

VCF_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie/Imputation/minDP3
cd /work-zfs/abattle4/heyuan/Variant_calling/GEO/QTL_calling
chromosome="$1"
	echo "chr${chromosome}"
	python construct_VCF.py ${chromosome}

	VCFfn=${VCF_dir}/dosage_by_sample_matrix_chr${chromosome}.recode.forFastQTL.vcf
	bgzip -c  ${VCFfn} > ${VCFfn}.gz
	tabix ${VCFfn}.gz
