#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=24:00:00
#SBATCH -p lrgmem

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis//ATAC_seq/alignment_bowtie

ml fastqtl
ml bcftools
ml htslib
ml python/3.7

# merge VCF files to get the correct header
cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie/VCF_files
merged_files=""
for f in `ls *forimputation.vcf.gz`
do
	echo $f
        tabix -f $f
        merged_files="${merged_files} ${f}"
done

bcftools merge ${merged_files} -o merged_vcf.vcf.gz -O z


tabix merged_vcf.vcf.gz
zcat merged_vcf.vcf.gz | grep "#" > merged_vcf_header.vcf
for chromosome in {1..22}
do
        cat merged_vcf_header.vcf > merged_vcf_chromosome${chromosome}.vcf
        tabix merged_vcf.vcf.gz ${chromosome} >> merged_vcf_chromosome${chromosome}.vcf
done

VCF_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/ATAC_seq/alignment_bowtie/Imputation/minDP3
cd /work-zfs/abattle4/heyuan/Variant_calling/GEO/QTL_calling
for chromosome in {1..22}
do
	echo "chr${chromosome}"
	python construct_VCF.py ${chromosome}

	VCFfn=${VCF_dir}/dosage_by_sample_matrix_chr${chromosome}.recode.forFastQTL.vcf
	bgzip -c  ${VCFfn} > ${VCFfn}.gz
	tabix ${VCFfn}.gz
done

