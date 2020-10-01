#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH -p shared

	for x in `ls /work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/ATAC_seq/alignment_bowtie/*.bowtie2.grch38.sortedByCoord.out.bam | xargs -n1 basename`
	do
		fn=${x/.bowtie2.grch38.sortedByCoord.out.bam/}
		echo $fn
		bash temp.sh ${fn}
	done
