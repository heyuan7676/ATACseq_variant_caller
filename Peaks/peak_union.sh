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
	cut -f1-3 ${f} | awk -vidx=$sample '{ print $0"\t"idx; }'  | sort -k1,1 -k2,2n > ${sample}.id.bed
done


# Take the union of all these ID-tagged files with BEDOPS bedops 
bedops --everything HG*.id.bed | bedmap --echo --echo-map-id-uniq --delim '\t' - > all.bed

# Get an interval that has regions in common with at least 5 input files 
awk -vthreshold=3 '(split($5, ids, ";") >= threshold)' all.bed | awk '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n | sed 's/ /	/g' > union-thresholded_level0.bed

# Merge intervals that overlap more than 50%
fn1=union-thresholded_level0.bed
for k in {1..10}
do
	echo level$k
	fn2=union-thresholded_level${k}.bed
	bedmap --count --echo-map-range --fraction-either 0.5 --delim '\t' ${fn1} \
    	| cut -f2- - \
    	| sort-bed - \
    	| uniq - \
    	> ${fn2}
	fn1=union-thresholded_level${k}.bed
done

wc -l union-thresholded_level*
rm -f *.count.unionPeaks.bed
rm -f *.count.unionPeaks.bed_matrix
rm -f peak_by_sample_matrix_chr*
rm -f union-peaks.bed
rm -f union-thresholded_level*clean*

