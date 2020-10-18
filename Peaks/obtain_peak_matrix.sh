#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6:00:00

cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks

header_fn=header.txt
echo "PeakID CHR START END" > ${header_fn}
for f in *count*bed
do
	#awk '{print $4,$5}' ${f} | sort -k1,1 > ${f}_temp
	sample=${f/.count.unionPeaks.bed/}
	echo ${sample} >> ${header_fn}
done

tr '\n' ' ' < ${header_fn} > ${header_fn}_temp
echo "" >>  ${header_fn}_temp


matrix_fn=peak_by_sample_matrix.txt
rm -f ${matrix_fn}
awk '{print $4, $1,$2,$3}' union-peaks.bed | sort -k1,1 >  ${matrix_fn}

for fn in `ls *count*bed_temp`
do
	echo $fn
	join -e0 -a 1 -a 2 -j 1 ${matrix_fn} -o auto ${fn} > temp
	mv temp ${matrix_fn}
done

sort -k2,2n -k3,3n ${matrix_fn} >> ${header_fn}_temp
mv ${header_fn}_temp ${matrix_fn}
awk '{print >> "peak_by_sample_matrix_chr"$2".txt"}' ${matrix_fn}  

rm *count*bed_temp
