#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6:00:00

cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Called_GT

## INFO
header_fn=header.txt
echo "CHR_POS" > ${header_fn}
for f in `ls ../VCF_files/*.filtered.recode.INFO.vcf | xargs -n1 basename`
do
	awk '{if($3>2) print $1"_"$2, $4}' ../VCF_files/${f} | sort -k1,1 -k2,2 > ${f}_temp
	sample=${f/.filtered.recode.INFO.vcf/}
	echo ${sample} >> ${header_fn}
done

tr '\n' ' ' < ${header_fn} > ${header_fn}_temp
echo "" >>  ${header_fn}_temp


unionFn=union-SNPs.bed
matrix_fn=gt_by_sample_matrix_INFO.txt
cat ${unionFn} > ${matrix_fn}

for fn in *.filtered.recode.INFO.vcf_temp
do
	echo $fn
	join -e0 -a 1 -a 2 -j 1 ${matrix_fn} -o auto ${fn} > temp
	mv temp ${matrix_fn}
done

sort -k2,2n -k3,3n ${matrix_fn} >> ${header_fn}_temp
mv ${header_fn}_temp ${matrix_fn}

rm *temp


