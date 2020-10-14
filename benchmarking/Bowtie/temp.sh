#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00

source ~/.bash_profile

cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks

# Count the number of reads overlapping with the peaks
ml bedtools
union_peaks=union-peaks.bed
for f in `ls *_peaks.narrowPeak | grep -v HG001`
do
        sample=${f/_peaks.narrowPeak/}
	echo $sample
	bam_file=../first_pass_bqsr/${sample}-clean.bam
	bedtools intersect -abam ${bam_file} -b ${union_peaks} -wo -bed | cut -d'	' -f13-16 | sort | uniq -c | awk '{print $2,$3,$4,$5,$1}' | sort -k1,1n -k2,2n | sed 's/ /	/g' > ${sample}.count.unionPeaks.bed
done

