#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH -p skylake

# Very fast and uses little memory

ml samtools
ml bcftools

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/
alignment="$1"


### STAR alignment
if [ "$alignment" == STAR ]
then
	alignment_dir=alignment_STAR
	variant_dir=${root_dir}/ATAC_seq/${alignment_dir}_variants
fi

### bowtie alignment
if [ "$alignment" == bowtie2 ]
then
	alignment_dir=alignment_bowtie
        variant_dir=${root_dir}/ATAC_seq/${alignment_dir}_variants
fi

ml vcftools
cd ${variant_dir}

# merge the vcf files and call genotype
rm -f merged*
str=""
for f in `ls *vcf.gz`; do str="${str} ${f}"; done
echo $str
bcftools merge $str -o merged.vcf
vcf-to-tab < merged.vcf > merged.raw.genotypes.tsv

