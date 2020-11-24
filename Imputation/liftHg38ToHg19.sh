#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH -p skylake

sample="$1"
cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/
VCF=${sample}.filtered.recode.vcf.gz

outputFn=GRCh37/${sample}.filtered.recode.GRCh37.vcf
gunzip ${VCF}


VCF=${VCF/.gz/}
cat ${VCF} |  grep "^#" > ${VCF}_temp
cat ${VCF} |  grep -v "^#" | awk '{print "chr"$0}' >> ${VCF}_temp

ml picard
picard LiftoverVcf I=${VCF}_temp O=${outputFn} CHAIN=/work-zfs/abattle4/heyuan/tools/liftOver/hg38ToHg19.over.chain.gz REJECT=${outputFn/.vcf/.rejected.vcf} R=/work-zfs/abattle4/heyuan/database/GRCh37_reference/hg19.fa
rm ${VCF}_temp


ml htslib
bgzip ${VCF}
