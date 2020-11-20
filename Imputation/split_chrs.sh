#!/bin/bash

sample="$1"
cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/
VCF=${sample}.filtered.recode.vcf.gz


ml htslib
ml vcftools
#bgzip -c $VCF #compress vcf
tabix -p vcf ${VCF} # index compressed vcf


mkdir -p ${sample}
tabix -H ${VCF} > ${sample}/${VCF}_header.txt
rm -f ${sample}/${sample}_chr*
for chromosome in {1..22}
do
        outfn=${sample}/${sample}_chr${chromosome}.vcf
        cat ${sample}/${VCF}_header.txt > ${outfn}
        tabix ${VCF} ${chromosome} | awk '{print "chr"$0}' >> ${outfn}
        bgzip ${outfn}
done

cd ${sample}
for f in *gz
do
	tabix -p vcf ${f}
done
