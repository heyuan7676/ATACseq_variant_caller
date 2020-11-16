#!/bin/bash

ml vcftools

VCF_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype

BOWTIE_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie
PEAK_DIR=${BOWTIE_DIR}/Peaks
peak_union_bed=${PEAK_DIR}/union-peaks.bed

OUT_DIR=${VCF_DIR}/WithInPeaks_Bowtie
mkdir -p ${OUT_DIR}

for chromosome in {1..22}
do
	vcf_file=${VCF_DIR}/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.GBRsampleKept.maf005.recode.vcf
	output_fn=${OUT_DIR}/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.GBRsampleKept.maf005.WithInPeaks

	vcftools --vcf ${vcf_file} --bed ${peak_union_bed} --out ${output_fn} --recode --keep-INFO-all

	# obtain real genotype instead of 0/1
        vcf-to-tab < ${output_fn}.recode.vcf > ${output_fn}.genotypes.tsv
done

