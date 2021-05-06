#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=12:00:00
#SBATCH -p lrgmem


qtl_filename=window1000000.fastq.permutation.results.BH
qtl_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/fastQTL
qtl_file=${qtl_dir}/${qtl_filename}.txt
random_file=${qtl_dir}/${qtl_filename}.randomMatched.txt

paste <(cat ${random_file} | cut -d'_' -f 1)  <(cat ${random_file} | awk '{print $1}' | cut -d'_' -f 2) <(cat ${random_file}) | awk '{print "chr"$1,$2-1,$2+1,$0}' | sed "1d" | sort -k1,1 -k2,2n | sed 's/ /	/g ' > ${random_file}_sorted
paste <(cat ${qtl_file} | awk '{print $2}' | cut -d'_' -f 1)  <(cat ${qtl_file} | awk '{print $2}' | cut -d'_' -f 2) <(cat ${qtl_file}) | awk '{print "chr"$1,$2-1,$2+1,$0}' | sed "1d" | sort -k1,1 -k2,2n | sed 's/ /	/g ' > ${qtl_file}_sorted


## CADD
cadd_file=/work-zfs/abattle4/lab_data/genomic_annotation_data/1000Genome/1000G.hg38.tsv
mkdir -p ${qtl_dir}/cadd
awk 'NR==FNR{inFileA[$2]; next} ($1 in inFileA)' ${qtl_file} ${cadd_file} > ${qtl_dir}/cadd/${qtl_filename}.cadd.txt
awk 'NR==FNR{inFileA[$2]; next} ($1 in inFileA)' ${random_file} ${cadd_file} > ${qtl_dir}/cadd/${qtl_filename}.randomMatched.cadd.txt


## TFBS
tfbs_file=/work-zfs/abattle4/lab_data/tfbs/remap2018_all_macs2_hg38_v1_2.bed
mkdir -p ${qtl_dir}/tfbs

bedtools intersect  -a ${qtl_file}_sorted -b ${tfbs_file} -wb >  ${qtl_dir}/tfbs/${qtl_filename}.tfbs.txt
bedtools intersect  -a ${random_file}_sorted -b ${tfbs_file} -wb >  ${qtl_dir}/tfbs/${qtl_filename}.randomMatched.tfbs.txt


## ENCODE
dnase_enhancer=/work-zfs/abattle4/lab_data/roadmap/dnase_enhancer/all_peaks.tsv.gz
bedtools intersect  -a ${qtl_file}_sorted -b ${dnase_enhancer} -wb >  ${qtl_dir}/tfbs/${qtl_filename}.dnase_enhancer.txt
bedtools intersect  -a ${random_file}_sorted -b ${dnase_enhancer} -wb >  ${qtl_dir}/tfbs/${qtl_filename}.randomMatched.dnase_enhancer.txt


## promoter
dnase_promoter=/work-zfs/abattle4/lab_data/roadmap/dnase_promoter/all_peaks.tsv.gz
bedtools intersect  -a ${qtl_file}_sorted -b ${dnase_promoter} -wb >  ${qtl_dir}/tfbs/${qtl_filename}.dnase_promoter.txt
bedtools intersect  -a ${random_file}_sorted -b ${dnase_promoter} -wb >  ${qtl_dir}/tfbs/${qtl_filename}.randomMatched.dnase_promoter.txt
