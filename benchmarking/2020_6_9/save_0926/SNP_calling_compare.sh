#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=10:00:00
#SBATCH -p shared

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/2019_10_2020_6/GM23248/

#cd ${root_dir}/genotyping
#bgzip ENCFF498RTM.vcf
#tabix -p vcf ENCFF498RTM.vcf.gz

#cd ${root_dir}/DNase/alignment
#for fn in `ls *bcf`
#do
#	bgzip ${fn}
#	tabix -p vcf ${fn}.gz
#done


### get the intersection of called variants
#ml bcftools
#cd ${root_dir}/DNase/alignment
#for fn in `ls *bcf.gz`
#do
#	mkdir ${fn}_intersect
#	bcftools isec -p ${fn}_intersect ${fn} ${root_dir}/genotyping/ENCFF498RTM.vcf.gz
#	paste <(cat ${fn}_intersect/0002.vcf | grep -v "^##" | awk '{print $1,$2,$3,$4,$5,$10}' | cut -d':' -f1 | sed 's/#//g') <(cat ${fn}_intersect/0002.vcf | grep -v "##" | awk '{print $8}' | cut -d';' -f1 | cut -d'=' -f2 ) > ${fn%.gz}_overlap_SNPs
#done

cd ${root_dir}/DNase/alignment
for exp in ENCLB544CVI ENCLB358EUZ
do
	bcftools isec -p ${exp}_replicate_intersect ${exp}_L1.raw.bcf.gz ${exp}_L2.raw.bcf.gz
	paste <(cat ${exp}_replicate_intersect/0002.vcf | grep -v "^##" | awk '{print $1,$2,$3,$4,$5,$10}' | cut -d':' -f1 | sed 's/#//g') <(cat ${exp}_replicate_intersect/0002.vcf | grep -v "##" | awk '{print $8}' | cut -d';' -f1 | cut -d'=' -f2)  | sed 's/ /	/g' > ${exp}_replicate_overlap_SNPs
done



### formatting
#cd /work-zfs/abattle4/heyuan/ASI/GM23248/DNase/alignment
#for fn in `ls *bcf.gz`
#do
#	paste <(zcat ${fn} | grep -v "^##" | awk '{print $1,$2,$3,$4,$5,$10}' | cut -d':' -f1 | sed 's/#//g' ) <(zcat ${fn}| grep -v "##" | awk '{print $8}' | cut -d';' -f1 | cut -d'=' -f2)> ${fn%.gz}_SNPs
#done


