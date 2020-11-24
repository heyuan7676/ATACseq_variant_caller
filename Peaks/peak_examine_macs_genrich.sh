#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00

macs2_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks
genrich_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Genrich
output_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Examine

mkdir -p ${output_dir}
ml bedtools


for sample in `cat test_all.txt`
do
	echo $sample
	fn1=${macs2_dir}/${sample}_peaks.narrowPeak
	fn2=${genrich_dir}/${sample}_peaks.Genrich.narrowPeak
	bedtools intersect -a ${fn2} -b ${fn1} > ${output_dir}/${sample}_genrich_in_macs2.bed
	bedtools intersect -a ${fn1} -b ${fn2} > ${output_dir}/${sample}_macs2_in_genrich.bed

	# for counting basepairs
	sort -k1,1n -k2,2n ${fn1} | bedtools merge > ${output_dir}/${sample}_macs2_peak_regions.bed
	sort -k1,1n -k2,2n ${fn2} | bedtools merge > ${output_dir}/${sample}_genrich_peak_regions.bed

 	bedtools intersect -a ${fn1} -b ${fn2} | sort -k1,1n -k2,2n | bedtools merge > ${output_dir}/${sample}_macs2_genrich_overlap.bed
done

