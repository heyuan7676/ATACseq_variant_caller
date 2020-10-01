#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH -p skylake

# Very fast and uses little memory

ml samtools
ml bcftools

fn="$1"
root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/

alignment=bowtie2


### STAR alignment
if [ "$alignment" == STAR ]
then
	# chromosome: 1,2,3, ... X
	genome_fa_dir=/work-zfs/abattle4/heyuan/database/STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa

	alignment_dir=alignment_STAR
	bam_dir=${root_dir}/ATAC_seq/${alignment_dir}
	bamFile=${fn}Aligned.sortedByCoord.out.passWF.removedDup.bam

	variant_dir=${root_dir}/ATAC_seq/${alignment_dir}_variants
	mkdir -p ${variant_dir}
	variant_fn=${variant_dir}/${fn}_bcftools
fi


### bowtie alignment
if [ "$alignment" == bowtie2 ]
then
	# chromosome: chr1,chr2,chr3, ... chgX
        genome_fa_dir=/work-zfs/abattle4/heyuan/database/GRCh38_reference/GRCh38.primary_assembly.genome.fa

        alignment_dir=alignment_bowtie
        bam_dir=${root_dir}/ATAC_seq/${alignment_dir}
        bamFile=${fn}.bowtie2.grch38.sortedByCoord.out.bam

        variant_dir=${root_dir}/ATAC_seq/${alignment_dir}_variants
        mkdir -p ${variant_dir}
        variant_fn=${variant_dir}/${fn}_bcftools
fi

# obtain real genotype instead of 0/1
ml vcftools
cd ${variant_dir}

ml htslib
bgzip ${variant_fn}.raw.vcf
bcftools index ${variant_fn}.raw.vcf.gz


# merge the vcf files and call genotype
#rm -f merged*
#str=""
#for f in `ls *vcf.gz`; do str="${str} ${f}"; done
#echo $str
#bcftools merge $str -o merged.vcf
#vcf-to-tab < merged.vcf > merged.raw.genotypes.tsv

