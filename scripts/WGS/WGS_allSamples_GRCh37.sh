#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH -p lrgmem

### Obtain the VCF file for the tested samples
onek_genome_dir=/work-zfs/abattle4/lab_data/1k_genomes

# cd ${onek_genome_dir}
# bash merge_VCF.sh

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/
mkdir -p ${root_dir}/Genotype
ml vcftools
ml bcftools

chromosome="$1"
echo chr$chromosome
fn=${onek_genome_dir}/ALL.chr${chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

keepSamples_fn=ALL.chr${chromosome}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.GBRsampleKept.maf005
vcftools \
	--gzvcf ${fn} \
	--keep samples_GBR.txt \
	--recode --recode-INFO-all \
	--maf 0.05 \
	--out ${root_dir}/Genotype/${keepSamples_fn}

cd ${root_dir}/Genotype
# obtain real genotype instead of 0/1
vcf-to-tab < ${keepSamples_fn}.recode.vcf > 1k_genome_chr${chromosome}.GRCh37.genotypes.tsv

# obtain MAF
bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\n' ${keepSamples_fn}.recode.vcf > 1k_genome_chr${chromosome}.GRCh37.AFs.tab

#gzip ${keepSamples_fn}.recode.vcf

