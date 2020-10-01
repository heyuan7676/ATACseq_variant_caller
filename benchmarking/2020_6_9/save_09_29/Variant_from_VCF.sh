#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH -p skylake

### Obtain the VCF file for the tested samples
onek_genome_dir=/work-zfs/abattle4/lab_data/1k_genomes_GRCh38/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL
# cd ${onek_genome_dir}
# bash merge_VCF.sh

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/
mkdir -p ${root_dir}/Genotype
ml vcftools
ml bcftools

chromosome="$1"
        echo chr$chromosome
        fn=${onek_genome_dir}/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

        keepSamples_fn=ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.GBRsampleKept.maf005
        vcftools \
        --gzvcf ${fn} \
        --keep samples_GBR.txt \
        --recode --recode-INFO-all \
        --maf 0.05 \
        --out ${root_dir}/Genotype/${keepSamples_fn}

        cd ${root_dir}/Genotype
        # obtain real genotype instead of 0/1
        vcf-to-tab < ${keepSamples_fn}.recode.vcf > 1k_genome_chr${chromosome}.genotypes.tsv

        # obtain MAF
        bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\n' ${keepSamples_fn}.recode.vcf > 1k_genome_chr${chromosome}.AFs.tab

	#gzip ${keepSamples_fn}.recode.vcf

