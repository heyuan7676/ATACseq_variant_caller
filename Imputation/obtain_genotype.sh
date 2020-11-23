#!/bin/bash

sample="$1"
datadir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/${sample}
liftOverdir=/work-zfs/abattle4/heyuan/tools/liftOver

ml vcftools
for chromosome in {1..22}
do
	cd ${datadir}
	#vcf-to-tab < chr${chromosome}.dose.vcf > chr${chromosome}.genotype.txt
	awk '{print "chr"$1,$2,$3 = $2 + 1,$4}' chr${chromosome}.genotype.txt | sed "1d" > chr${chromosome}.genotype.bed

	cd ${liftOverdir}
	./liftOver ${datadir}/chr${chromosome}.genotype.bed hg19ToHg38.over.chain.gz ${datadir}/chr${chromosome}.genotype.GRCh38.txt unlifted.bed  
	
done

cd ${datadir}
head -n1 chr${chromosome}.genotype.txt > ${sample}.imputed.genotype.GRCh38.txt
cat chr*.genotype.GRCh38.txt | sed 's/chr//g' >> ${sample}.imputed.genotype.GRCh38.txt
