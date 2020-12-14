#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00


sample1="$1"
sample2="$2"
ROOT_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/PRJNA603385/ATAC_seq/alignment_bowtie/
macs2_dir=${ROOT_DIR}/Peaks
output_dir=${ROOT_DIR}/Peaks_Examine

fn1=${macs2_dir}/${sample1}_peaks.narrowPeak
fn2=${macs2_dir}/${sample2}_peaks.narrowPeak

mkdir -p ${output_dir}
ml bedtools
bedtools intersect -a ${fn2} -b ${fn1} > ${output_dir}/${sample1}_${sample2}_macs2_overlap.bed

# for counting basepairs
sort -k1,1n -k2,2n ${fn1} | bedtools merge > ${output_dir}/${sample1}_macs2_peak_regions.bed
sort -k1,1n -k2,2n ${fn2} | bedtools merge > ${output_dir}/${sample2}_macs2_peak_regions.bed
bedtools intersect -a ${fn1} -b ${fn2} | sort -k1,1n -k2,2n | bedtools merge > ${output_dir}/${sample1}_${sample2}_macs2_overlap_bp.bed
rm ${output_dir}/${sample1}_macs2_peak_regions.bed
rm ${output_dir}/${sample2}_macs2_peak_regions.bed

