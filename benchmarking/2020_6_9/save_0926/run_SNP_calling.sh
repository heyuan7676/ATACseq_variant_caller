#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH -p shared

alignment="$1"

### STAR alignment
if [ "$alignment" == STAR ]
then
	for x in `ls /work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/ATAC_seq/alignment_STAR/*removedDup.bam | xargs -n1 basename`
	do 
		fn=${x/Aligned.sortedByCoord.out.passWF.removedDup.bam/}
		echo $fn
		bash variantCaller_bcftools.sh ${fn}
	done
fi



### STAR alignment
if [ "$alignment" == bowtie2 ]
then
	for x in `ls /work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/ATAC_seq/alignment_bowtie/*.bowtie2.grch38.sortedByCoord.out.bam | xargs -n1 basename`
	do
		fn=${x/.bowtie2.grch38.sortedByCoord.out.bam/}
		echo $fn
		sbatch variantCaller_bcftools.sh ${fn} bowtie2
	done
fi
