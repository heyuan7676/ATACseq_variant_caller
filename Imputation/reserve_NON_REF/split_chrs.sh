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
for chromosome in {22..22}
do
        fn_chr=${sample}/${sample}_chr${chromosome}_temp.vcf
        tabix ${VCF} ${chromosome} | awk '{print "chr"$0}' > ${fn_chr}
	python /work-zfs/abattle4/heyuan/Variant_calling/Imputation/substitute_NON_REF.py ${fn_chr}
done

cd ${sample}
for chromosome in {22..22}
do
	outfn=${sample}_chr${chromosome}.vcf
	fn_chr=${sample}_chr${chromosome}_temp.vcf

	cat ${VCF}_header.txt > ${outfn}	
	sed "1d" ${fn_chr} | sed 's/,<NON_REF>//g' >> ${outfn}
	rm ${fn_chr}

	bgzip ${outfn}
	tabix -p vcf ${outfn}.gz
done
