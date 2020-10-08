#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6:00:00


ml python/3.7
variant_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/ATAC_seq/alignment_bowtie/
for sample in `ls ${variant_dir}*DP9.txt | xargs -n1 basename | cut -d'-' -f1`
do
	echo $sample
	python Evaluation_metrics.py ${sample}
done
