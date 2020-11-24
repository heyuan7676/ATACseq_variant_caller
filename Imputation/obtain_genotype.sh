#!/bin/bash

sample="$1"
datadir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/${sample}
liftOverdir=/work-zfs/abattle4/heyuan/tools/liftOver

ml vcftools
for chromosome in {1..22}
do
	cd ${datadir}
	#vcf-to-tab < chr${chromosome}.dose.vcf > chr${chromosome}.genotype.txt
	awk '{print "chr"$1,$2,$3 = $2 + 1,$4}' chr${chromosome}.imputed.GRCh37.genotype.txt | sed "1d" > chr${chromosome}.imputed.GRCh37.genotype.bed

	cd ${liftOverdir}
	./liftOver ${datadir}/chr${chromosome}.imputed.GRCh37.genotype.bed hg19ToHg38.over.chain.gz ${datadir}/chr${chromosome}.imputed.GRCh38.genotype.txt unlifted.bed  
	
done

cd ${datadir}
head -n1 chr${chromosome}.imputed.GRCh37.genotype.txt > ${sample}.imputed.GRCh38.genotype.txt
cat chr*.imputed.GRCh38.genotype.txt | sed 's/chr//g' >> ${sample}.imputed.GRCh38.genotype.txt
