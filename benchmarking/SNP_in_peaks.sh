#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH -p skylake

ml bedtools

exp_fn="$1"
#exp_fn=ENCLB358EUZ
#exp_fn=ENCLB544CVI

vcfFn=/work-zfs/abattle4/heyuan/ASI/GM23248/genotyping/ENCFF498RTM.vcf
pFn=/work-zfs/abattle4/heyuan/ASI/GM23248/DNase/peaks/${exp_fn}.peaks.bed
outdir=/work-zfs/abattle4/heyuan/ASI/GM23248/DNase/SNP_calling/

intersectBed -a ${vcfFn} -b ${pFn} >  ${outdir}/WGS_SNPs_in_${exp_fn}.vcf
cat ${outdir}/WGS_SNPs_in_${exp_fn}.vcf | grep "0/1" > ${outdir}/WGS_SNPs_in_${exp_fn}_heterozygous.vcf

### SNP from DNase-seq

ml samtools
ml python/2.7
alignmentdir=/work-zfs/abattle4/heyuan/ASI/GM23248/DNase/alignment
for fn in `ls ${alignmentdir}/${exp_fn}*bam | xargs -n1 basename`
do
	echo $fn
	intersectBed -a ${alignmentdir}/${fn} -b ${pFn} > ${outdir}/DNase_SNPs_in_${exp_fn}_${fn}
	samtools view ${outdir}/DNase_SNPs_in_${exp_fn}_${fn} | grep vA | grep vW:i:1 | grep -v "vA:B:c,3" | sed 's/CHR//g' > ${outdir}/DNase_SNPs_in_${exp_fn}_${fn}.filtered.txt

	python format_filtered_bam.py ${outdir}/DNase_SNPs_in_${exp_fn}_${fn}.filtered.txt
done

