#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=6:00:00
#SBATCH -p shared

STAR_DIR=/work-zfs/abattle4/heyuan/tools/STAR-2.7.1a/bin/Linux_x86_64_static

### 0. Prepare genome index
#cd ${STAR_DIR}
#genome_fa_dir=/work-zfs/abattle4/heyuan/database/STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa
#gtf_dir=/work-zfs/abattle4/heyuan/database/STAR/gencode_chromosomenames.v30.annotation.gtf
#genome_index_idr=/work-zfs/abattle4/heyuan/database/STAR/genome_index

#./STAR  --runMode genomeGenerate --runThreadN 6  --genomeDir ${genome_index_idr} --genomeFastaFiles ${genome_fa_dir} --sjdbGTFfile ${gtf_dir} --limitGenomeGenerateRAM 115999096192 

fn="$1"
#fn=HG00246
echo $fn

genome_index_idr=/work-zfs/abattle4/heyuan/database/STAR/genome_index

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/
#vcfFn=${root_dir}/vcf_Fn/ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.GBRsampleRemoved.maf005.recode.vcf
vcfFn=${root_dir}/vcf_Fn/ALL.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.GBRsampleRemoved.maf005.recode.vcf

fastq_dir=${root_dir}/ATAC_seq/ftp.sra.ebi.ac.uk/vol1/ERZ683/
fn1=${fn}_r1.fixed.fastq.gz
fn2=${fn}_r2.fixed.fastq.gz

alignment_dir=${root_dir}/ATAC_seq/alignment_STAR/
mkdir -p ${alignment_dir}
alignment_fn=${alignment_dir}/${fn}

final_bamFile=${alignment_dir}/${fn}Aligned.sortedByCoord.out.passWF.removedDup.bam
if [ -f ${final_bamFile} ]
then
	echo "${final_bamFile} exists"
	exit
fi

### 1. Align reads to the genome
cd ${STAR_DIR}
./STAR --genomeDir ${genome_index_idr}  \
	--runThreadN 12 \
	--readFilesIn ${fastq_dir}/${fn1} ${fastq_dir}/${fn2} \
	--readFilesCommand zcat  \
	--varVCFfile ${vcfFn}  \
	--outFileNamePrefix ${alignment_fn} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMattributes NH HI AS nM NM MD vA vG vW \
	--waspOutputMode SAMtag \
	--bamRemoveDuplicatesType UniqueIdentical \
	--limitBAMsortRAM 35690842466

# output
bamFile=${alignment_fn}Aligned.sortedByCoord.out.bam

### 2. Preserve reads that pass the WASP filter (automatically include only reads that map to variants)
# vW: WASP filtering for unbiased allele-specific read mapping <-- especially important in this task
ml samtools
output_bamFile=${alignment_dir}/${fn}Aligned.sortedByCoord.out.passWF.bam
samtools view -H $bamFile > ${output_bamFile}_header
samtools view ${bamFile} | grep "vW:i:1" | cat ${output_bamFile}_header - | samtools view -b -o ${output_bamFile}
rm ${output_bamFile}_header
rm ${bamFile}

## Remove duplicate reads
final_bamFile=${alignment_dir}/${fn}Aligned.sortedByCoord.out.passWF.removedDup.bam
samtools view -b -F 0x400 ${output_bamFile} > ${final_bamFile}
rm ${output_bamFile}




