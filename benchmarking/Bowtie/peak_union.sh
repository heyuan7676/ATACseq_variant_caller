#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6:00:00

source ~/.bash_profile

cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks

# Label your input BED files so that their IDs uniquely identify their intervals
for f in *_peaks.narrowPeak
do
	sample=${f/_peaks.narrowPeak/}
	cut -f1-3 ${f} | awk -vidx=$sample '{ print $0"\t"idx; }' > ${sample}.id.bed
done


# Take the union of all these ID-tagged files with BEDOPS bedops 
bedops --everything *.id.bed | bedmap --echo --echo-map-id-uniq --delim '\t' - > all.bed

# Get an interval that has regions in common with at least 5 input files 
awk -vthreshold=3 '(split($5, ids, ";") >= threshold)' all.bed | awk '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n | sed 's/ /	/g' > union-thresholded_level0.bed

# Merge intervals that overlap more than 50%
fn1=union-thresholded_level0.bed
for k in {1..10}
do
	echo level$k
	fn2=union-thresholded_level${k}.bed
	bedmap --count --echo-map-range --fraction-both 0.5 --delim '\t' ${fn1} \
    	| cut -f2- - \
    	| sort-bed - \
    	| uniq - \
    	> ${fn2}
	fn1=union-thresholded_level${k}.bed
done

wc -l union-thresholded_level*

# Use union-thresholded_level10.bed to count number of reads in each peak
union_peaks=union-peaks.bed
awk '{print $0,"Peak"NR}' union-thresholded_level${k}.bed | sed 's/ /	/g' > ${union_peaks}

# Count the number of reads overlapping with the peaks
ml bedtools
for f in *_peaks.narrowPeak
do
        sample=${f/_peaks.narrowPeak/}
	echo $sample
	bam_file=../first_pass_bqsr/${sample}-clean.bam
	bedtools intersect -abam ${bam_file} -b ${union_peaks} -wo -bed | cut -d'	' -f13-16 | sort | uniq -c | awk '{print $2,$3,$4,$5,$1}' | sort -k1,1n -k2,2n | sed 's/ /	/g' > ${sample}.count.unionPeaks.bed
done

