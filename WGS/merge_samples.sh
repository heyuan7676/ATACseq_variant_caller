#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH -p skylake

ml vcftools
ml bcftools
ml samtools
ml htslib

VCF_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype

BOWTIE_DIR=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie
PEAK_DIR=${BOWTIE_DIR}/Peaks


cd ${VCF_DIR}/WithInPeaks_Bowtie
for chromosome in {1..22}
do
        rm -f chr${chromosome}.WithInPeaks.vcf
        rm -f chr${chromosome}.WithInPeaks.genotypes.tsv

        str=""
        for f in `ls HG*/*.chr${chromosome}.WithInPeaks.recode.vcf.gz`; do
                str="${str} ${f}"
        done
        echo $str

        bcftools merge $str -o chr${chromosome}.WithInPeaks.vcf

        vcf-to-tab < chr${chromosome}.WithInPeaks.vcf > chr${chromosome}.WithInPeaks.genotypes.tsv
done
