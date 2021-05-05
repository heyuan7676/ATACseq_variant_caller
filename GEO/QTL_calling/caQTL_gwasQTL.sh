#!/bin/bash
#SBATCH --time 50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --no-requeue

### produce:  SNPs x Tissues in DNase data



gwas_dir=/work-zfs/abattle4/parsana/gtex_trans_v8/data/gwas_haky/significant_5e8
caQTL_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/fastQTL
outdir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/Haematopoiesis/Gencove/coloc_gwasQTLs

PREFIX=global

mkdir -p ${outdir}
caQTL_file=window1000000.fastq.permutation.results.BH.txt
awk '{if($7<0.05) print $2"_"}' ${caQTL_dir}/${caQTL_file} > ${outdir}/${caQTL_file/.txt/.variants.txt}

cd ${gwas_dir}
for fn in `ls *txt`
do
	N=`wc -l ${fn} | awk '{print $1}'`
	if [ "$N" -lt 100 ]
	then
		continue
	fi
	echo $fn
	outfn=${PREFIX}_gwas_${fn}
	grep -Ff ${outdir}/${caQTL_file/.txt/.variants.txt} ${fn} | awk '{print $2,$0}'> ${outdir}/${outfn}
done	



eQTL_dir=/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL

cd ${eQTL_dir}
for fn in `ls *.v8.signif_variant_gene_pairs.txt`
do
	echo $fn
	outfn=${PREFIX}_eQTL_${fn}
	grep -Ff ${outdir}/${caQTL_file/.txt/.variants.txt} ${fn} > ${outdir}/${outfn}
done	


cd ${outdir}
for file1 in ${PREFIX}*_gwas_*
do
	for file2 in ${PREFIX}*_eQTL_*
	do
		join <(sort $file1) <(sort $file2) > test_middle.txt
		N=`wc -l test_middle.txt | awk '{print $1}'`
        	if [ "$N" -gt 0 ]
        	then
			mv test_middle.txt ${PREFIX}_both_${file1/.txt/}_${file2/.v8.signif_variant_gene_pairs.txt/}.txt
		fi
	done
done
