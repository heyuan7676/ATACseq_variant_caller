#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=24:00:00
#SBATCH -p shared

ml samtools
ml bcftools

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/

alignment="$1"


### STAR alignment
if [ "$alignment" == STAR ]
then
	# chromosome: 1,2,3, ... X
	genome_fa_dir=/work-zfs/abattle4/heyuan/database/STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa

	alignment_dir=alignment_STAR
	bam_dir=${root_dir}/ATAC_seq/${alignment_dir}

	cd ${bam_dir}
	str=""
	for f in `ls *Aligned.sortedByCoord.out.passWF.removedDup.bam`
	do
		if [[ ! -e ${1}.bai ]];
  		then
  			samtools index $f
		fi
		str="${str} ${f}"
	done
	bamFile=${str}
	echo $bamFile

	variant_dir=${root_dir}/ATAC_seq/${alignment_dir}_variants
	mkdir -p ${variant_dir}
	variant_fn=${variant_dir}/multipleSamples_bcftools
fi


### bowtie alignment
if [ "$alignment" == bowtie2 ]
then
	# chromosome: chr1,chr2,chr3, ... chgX
        genome_fa_dir=/work-zfs/abattle4/heyuan/database/GRCh38_reference/GRCh38.primary_assembly.genome.fa

        alignment_dir=alignment_bowtie
        bam_dir=${root_dir}/ATAC_seq/${alignment_dir}

	cd ${bam_dir}
        str=""
        for f in `ls *.bowtie2.grch38.sortedByCoord.out.bam`
        do
                if [[ ! -e ${1}.bai ]];
                then    
                        samtools index $f
                fi
                str="${str} ${f}"
        done
        bamFile=${str}
        echo $bamFile

        variant_dir=${root_dir}/ATAC_seq/${alignment_dir}_variants
        mkdir -p ${variant_dir}
        variant_fn=${variant_dir}/multipleSamples_bcftools
fi


# Q: Minimum base quality for a base to be considered
# A: Do not skip anomalous read pairs in variant calling
# Ou: output uncompressed BCF
# f: faidx-indexed reference file 

bcftools mpileup \
  -Q 30 \
  -A \
  -Ou \
  -f ${genome_fa_dir} \
  ${bamFile} | bcftools call -v -c -Ov - > ${variant_fn}.raw.vcf

# obtain real genotype instead of 0/1
ml vcftools
cd ${variant_dir}

ml htslib
bgzip ${variant_fn}.raw.vcf
bcftools index ${variant_fn}.raw.vcf.gz


