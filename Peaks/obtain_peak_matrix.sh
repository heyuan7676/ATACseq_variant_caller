#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6:00:00

cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks

for f in *count*bed
do
	awk '{print $4,$5}' ${f} | sort -k1,1 > ${f}_temp
done


matrix_fn=peak_by_sample_matrix
cat ${f}_temp > ${matrix_fn}

for f in *count*bed_temp
do
	join -a 1 -a 2 -j 1 ${f} ${matrix_fn} > temp
	mv temp ${matrix_fn}
done

rm *count*bed_temp
