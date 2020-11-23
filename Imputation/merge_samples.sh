#!/bin/bash


cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files

ml bcftools

for chromosome in {22..22}
do
	str=""
	for t in `cat test_all.txt | head -n2`
	do 
		str=${str}" "${t}/${t}_chr${chromosome}.vcf.gz
	done
	outfn=chr${chromosome}_allsamples.forMIS.vcf
	bcftools merge ${str} > ${outfn}
	bgzip ${outfn}
done
