#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00

ROOT_DIR="$1"
genrich_dir=${ROOT_DIR}/Peaks_Genrich
combined_dir=${ROOT_DIR}/Peaks_Genrich/combined
output_dir=${ROOT_DIR}/Peaks_Examine

fn1=${genrich_dir}/union-peaks.bed
fn2=${combined_dir}/union-peaks_combined.bed 

mkdir -p ${output_dir}
ml bedtools
bedtools intersect -a ${fn2} -b ${fn1} > ${output_dir}/combined_in_genrich.bed
bedtools intersect -a ${fn1} -b ${fn2} > ${output_dir}/genrich_in_combined.bed

# for counting basepairs
sort -k1,1n -k2,2n ${fn1} | bedtools merge > ${output_dir}/genrich_peak_regions.bed
sort -k1,1n -k2,2n ${fn2} | bedtools merge > ${output_dir}/combined_peak_regions.bed
bedtools intersect -a ${fn1} -b ${fn2} | sort -k1,1n -k2,2n | bedtools merge > ${output_dir}/genrich_combined_overlap.bed

