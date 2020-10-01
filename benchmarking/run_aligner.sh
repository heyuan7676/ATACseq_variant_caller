#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=12:00:00
#SBATCH -p shared


for f in `ls /work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/ATAC_seq/ftp.sra.ebi.ac.uk/vol1/ERZ683/*r2*gz | xargs -n1 basename`
do 
	sampleID=${f/_r2.fixed.fastq.gz/}
	#fn=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR//ATAC_seq/alignment_bowtie/${sampleID}.bowtie2.grch38.sortedByCoord.out.bam
	#if [ -f ${fn} ]
	#then
#		echo $sampleID
#		sbatch aligner_bowtie2.sh ${sampleID}
#	else
#		continue
#	fi
	sbatch aligner_bowtie2.sh ${sampleID}
	#bash aligner_STAR.sh ${f}
done
