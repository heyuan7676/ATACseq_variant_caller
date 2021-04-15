#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH -p lrgmem

peak_calling=MACS2
method=Imputation

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis//ATAC_seq/alignment_bowtie

ml fastqtl
ml bcftools
ml htslib


chromosome="$1"
#for chromosome in {1..22}
#do
	peakDir=${root_dir}/Peaks_${peak_calling}
	peakFn=${peakDir}/peak_by_sample_matrix_RPKM_corrected_chromosome${chromosome}_forFastQTL.bed

	rm -f  ${peakFn}.gz
	bgzip ${peakFn}&& tabix -p bed ${peakFn}.gz

	VCF_dir=${root_dir}/${method}/minDP3
	VCFfn=${VCF_dir}/dosage_by_sample_matrix_chr${chromosome}.recode.forFastQTL.vcf

	outdir=${root_dir}/fastQTL/${method}/${peak_calling}
	mkdir -p ${outdir}

	#for cisDis in 500 1000 10000
	for cisDis in 1000 10000
	do
		fastqtl --vcf ${VCFfn}.gz --bed ${peakFn}.gz --window ${cisDis} --region chr${chromosome}:1-300000000 --out ${outdir}/chr${chromosome}.window${cisDis}.fastq.results
		#fastqtl --vcf ${VCFfn}.gz --bed ${peakFn}.gz --window ${cisDis}  --region chr${chromosome}:1-300000000 --permute 1000 --out ${outdir}/chr${chromosome}.window${cisDis}.fastq.permutation.results
	done

#done
