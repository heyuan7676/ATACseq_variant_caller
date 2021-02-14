#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH -p skylake

peak_calling="$1"
chromosome="$2"

outdir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/fastQTL/${peak_calling}
mkdir -p ${outdir}
ml fastqtl


VCFfn=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/ALL.chr${chromosome}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.vcf
#cat ${VCFfn} | grep "#" > ${VCFfn}_header
#cat ${VCFfn} | grep -v "#" | awk '{$3=$1"_"$2; print "chr"$0}' >> ${VCFfn}_header 
#cat ${VCFfn}_header | sed 's/ /	/g' > ${VCFfn}.withchr

#ml htslib
#bgzip -c  ${VCFfn}.withchr > ${VCFfn}.withchr.gz
#tabix ${VCFfn}.withchr.gz


ml htslib
peakDir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_${peak_calling}
peakFn=${peakDir}/peak_by_sample_matrix_RPKM_corrected_chromosome${chromosome}_forFastQTL.bed

#rm -f  ${peakFn}.gz
bgzip ${peakFn}&& tabix -p bed ${peakFn}.gz


for cisDis in 1 500 1000 10000 100000 1000000
do
    fastqtl --vcf ${VCFfn}.withchr.gz --bed ${peakFn}.gz --window ${cisDis} --region chr${chromosome}:1-300000000 --out ${outdir}/chr${chromosome}.window${cisDis}.fastq.results
    fastqtl --vcf ${VCFfn}.withchr.gz --bed ${peakFn}.gz --window ${cisDis}  --region chr${chromosome}:1-300000000 --permute 1000 --out ${outdir}/chr${chromosome}.window${cisDis}.fastq.permutation.results
done



