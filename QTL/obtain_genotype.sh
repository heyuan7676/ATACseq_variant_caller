#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6:00:00

cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Called_GT

header_fn=header.txt
echo "CHR POS" > ${header_fn}
for f in *filtered.genotype.minDP3.txt
do
	awk '{print $1"_"$2,$4}' ${f} | sed 1d | sort -k1,1 > ${f}_temp
	sample=${f/.filtered.genotype.minDP3.txt/}
	echo ${sample} >> ${header_fn}
done

tr '\n' ' ' < ${header_fn} > ${header_fn}_temp
echo "" >>  ${header_fn}_temp


unionFn=union-SNPs.bed
rm -f {unionFn}
for f in *filtered.genotype.minDP3.txt_temp
do
	echo $f
	awk '{print $1}' ${f} >> ${unionFn}
	sort ${unionFn} | uniq > ${unionFn}_temp
	mv ${unionFn}_temp ${unionFn}
done

matrix_fn=gt_by_sample_matrix.txt
cat ${unionFn} > ${matrix_fn}

for fn in *filtered.genotype.minDP3.txt_temp
do
	echo $fn
	join -e0 -a 1 -a 2 -j 1 ${matrix_fn} -o auto ${fn} > temp
	mv temp ${matrix_fn}
done

sort -k2,2n -k3,3n ${matrix_fn} >> ${header_fn}_temp
mv ${header_fn}_temp ${matrix_fn}

rm *count*bed_temp
