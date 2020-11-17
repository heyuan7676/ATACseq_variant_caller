#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH -p skylake

ml vcftools
ml bcftools
ml samtools
ml htslib

VCF_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype

BOWTIE_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie
PEAK_DIR=${BOWTIE_DIR}/Peaks

sample="$1"
peak_file=${PEAK_DIR}/${sample}_peaks.narrowPeak

str=""
for chromosome in {1..22}
do
        vcf_file=${VCF_DIR}/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.GBRsampleKept.maf005.recode.vcf

	OUT_DIR=${VCF_DIR}/WithInPeaks_Bowtie/${sample}
	mkdir -p ${OUT_DIR}
	output_fn=${OUT_DIR}/${sample}.chr${chromosome}.WithInPeaks

	vcftools \
	--vcf ${vcf_file} \
	--bed ${peak_file} \
	--indv ${sample} \
	--out ${output_fn} \
	--recode \
	--keep-INFO-all

	output_fn=${OUT_DIR}/${sample}.chr${chromosome}.WithInPeaks.recode.vcf
	bgzip ${output_fn}
	bcftools index ${output_fn}.gz

	str="${str} ${output_fn}.gz"

done

echo $str
bcftools concat $str -o ${OUT_DIR}/${sample}.WithInPeaks.vcf
vcf-to-tab < ${OUT_DIR}/${sample}.WithInPeaks.vcf > ${OUT_DIR}/${sample}.WithInPeaks.genotypes.tsv


