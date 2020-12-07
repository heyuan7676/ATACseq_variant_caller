#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=6:00:00

### Obtain the VCF file for the tested samples
onek_genome_dir=/work-zfs/abattle4/lab_data/1k_genomes_GRCh38/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/

# cd ${onek_genome_dir}
# bash merge_VCF.sh

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/
mkdir -p ${root_dir}/Genotype
ml vcftools
ml bcftools

chromosome="$1"
echo chr$chromosome
fn=${onek_genome_dir}/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

keepMAF_fn=ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005
vcftools --gzvcf ${fn} \
	--maf 0.05 \
	--recode --recode-INFO-all \
	--out ${root_dir}/Genotype/${keepMAF_fn}


cd ${root_dir}/Genotype

echo "CHR_POS #CHROM POS" > ${keepMAF_fn}.variants.txt
awk '{print $1"_"$2, $1,$2}' ${keepMAF_fn}.recode.vcf | grep -v "#" >> ${keepMAF_fn}.variants.txt

rm ${keepMAF_fn}.recode.vcf
