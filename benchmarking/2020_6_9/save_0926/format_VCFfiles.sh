#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH -p skylake

### Obtain the VCF file excluding tested samples
onek_genome_dir=/work-zfs/abattle4/lab_data/1k_genomes_GRCh38/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL
# cd ${onek_genome_dir}
# bash merge_VCF.sh

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/
mkdir -p ${root_dir}/vcf_Fn
ml vcftools 

chromosome="$1"
	echo chr$chromosome
	fn=${onek_genome_dir}/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
	fn_v43=${onek_genome_dir}/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.v43.vcf.gz
        #mv ${fn} ${fn_v43}
        #echo $fn
        #zcat $fn_v43 | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > ${fn/.gz/}
        #gzip ${fn/.gz/}

	SampleRemoved_fn=ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.GBRsampleRemoved.maf005
	vcftools \
	--gzvcf ${fn} \
	--remove samples_GBR.txt \
	--recode --recode-INFO-all \
	--maf 0.05 \
	--out ${root_dir}/vcf_Fn/${SampleRemoved_fn}

	cd ${root_dir}/vcf_Fn
	gzip ${SampleRemoved_fn}.recode.vcf

