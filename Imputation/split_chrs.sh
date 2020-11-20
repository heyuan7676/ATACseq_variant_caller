#!/bin/bash

sample="$1"
cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/
VCF=${sample}.filtered.recode.vcf.gz


ml htslib
ml vcftools
ml python/2.7

#bgzip -c $VCF #compress vcf
tabix -p vcf ${VCF} # index compressed vcf


mkdir -p ${sample}
tabix -H ${VCF} > ${sample}/${VCF}_header.txt
rm -f ${sample}/${sample}_chr*
for chromosome in {1..22}
do
	outfn=${sample}/${sample}_chr${chromosome}.vcf
	cat ${sample}/${VCF}_header.txt > ${outfn}
        tabix ${VCF} ${chromosome} | awk '{print "chr"$0}' | sed '/	<NON_REF>/d' | sed 's/,<NON_REF>//g' >> ${outfn}
done

cd ${sample}
for chromosome in {1..22}
do
	outfn=${sample}_chr${chromosome}.vcf
	bgzip ${outfn}
        tabix -p vcf ${outfn}.gz
	bcftools annotate -x INFO ${outfn}.gz |  bcftools annotate -x FORMAT > temp
	rm ${outfn}.gz*
	mv temp ${outfn}
	bgzip ${outfn}
	tabix -p vcf ${outfn}.gz
done
