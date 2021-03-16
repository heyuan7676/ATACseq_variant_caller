#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6:00:00


ROOT_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie
cd ${ROOT_DIR}/Peaks_Genrich

ml bedtools

# Label your input BED files so that their IDs uniquely identify their intervals
for f in *_peaks.Genrich.narrowPeak 
do
	sample=${f/_peaks.Genrich.narrowPeak/}
	echo $sample
        bedtools intersect -wao -a union-peaks.bed -b ${f} > ${sample}_in_union_peaks.bed
done

