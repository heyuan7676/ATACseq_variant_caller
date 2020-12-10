#!/bin/bash
#SBATCH --time 20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --no-requeue
#SBATCH -p lrgmem

#bash vcfFn_process.sh

## compute LD blocks
module load plink
chromosome="$1"
r2Thr=0.2


onek_genome_dir=/work-zfs/abattle4/lab_data/1k_genomes_GRCh38/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/

# cd ${onek_genome_dir}
# bash merge_VCF.sh

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/
ml vcftools
ml bcftools

echo chr$chromosome
fn=${onek_genome_dir}/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
mkdir -p ${root_dir}/Genotype
mkdir -p ${root_dir}/plink

keepMAF_fn=ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005
vcftools --gzvcf ${fn} \
        --maf 0.05 \
        --recode --recode-INFO-all \
        --out ${root_dir}/Genotype/${keepMAF_fn}

## transform VCF files to plink binary files
vcffn=${root_dir}/Genotype/${keepMAF_fn}.recode.vcf
binaryFn=${root_dir}/plink/${keepMAF_fn}.plink
plink --vcf ${vcfFn_chr} --recode --out ${binaryFn}


## compute R2
r2Out=${root_dir}/plink/${keepMAF_fn}.plink.r2
plink --file ${binaryFn} --r2 inter-chr --ld-window-r2 ${r2Thr} --out ${r2Out}

