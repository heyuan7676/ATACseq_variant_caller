#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00

# Very fast and uses little memory

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

cd ${variant_dir}

for f in `ls *vcf.gz`
do
	echo $f
	bcftools query -f '%CHROM\t%POS\tDP=%DP\tAF1=%AF1\tAC1=%AC1\n' ${f} > ${f/.vcf.gz/.vcf.INFO.txt}
done

