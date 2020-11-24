#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00

macs2_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks
combined_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks/combined
output_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Examine


fn1=${macs2_dir}/union-peaks.bed
fn2=${macs2_dir}/combined/union-peaks_combined.bed 

mkdir -p ${output_dir}
ml bedtools
bedtools intersect -a ${fn2} -b ${fn1} > ${output_dir}/combined_in_macs2.bed
bedtools intersect -a ${fn1} -b ${fn2} > ${output_dir}/macs2_in_combined.bed

# for counting basepairs
sort -k1,1n -k2,2n ${fn1} | bedtools merge > ${output_dir}/macs2_peak_regions.bed
sort -k1,1n -k2,2n ${fn2} | bedtools merge > ${output_dir}/combined_peak_regions.bed
bedtools intersect -a ${fn1} -b ${fn2} | sort -k1,1n -k2,2n | bedtools merge > ${output_dir}/macs2_combined_overlap.bed

