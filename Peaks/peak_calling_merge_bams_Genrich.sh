#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=24:00:00
#SBATCH -p lrgmem


BOWTIE_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie

DIR_FIRST_PASS=${BOWTIE_DIR}/first_pass_bqsr
BLACKLIST=${BOWTIE_DIR}/Peaks/hg38.blacklist.bed

str=${DIR_FIRST_PASS}/`head -n1 test_all.txt`"-clean.sorted.bam"
for s in `cat test_all.txt | sed 1d`
do
	str=${DIR_FIRST_PASS}/${s}"-clean.sorted.bam,"${str}
done

echo $str

signal=${BOWTIE_DIR}/Peaks_Genrich/combined_peaks.Genrich.narrowPeak
bedgraph=${BOWTIE_DIR}/Peaks_Genrich/combined_peaks.Genrich.bedgraph
/work-zfs/abattle4/heyuan/tools/Genrich-master/Genrich -t ${str} -j -r -o ${signal} -f ${bedgraph} -v -E ${BLACKLIST}


