#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=10:00:00
#SBATCH -p skylake

ml bowtie2
ml samtools

genome_index=/work-zfs/abattle4/heyuan/database/bowtie2/grch38/grch38_1kgmaj

fn="$1"
echo $fn

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/
fastq_dir=${root_dir}/ATAC_seq/ftp.sra.ebi.ac.uk/vol1/ERZ683/
fn1=${fn}_r1.fixed.fastq.gz
fn2=${fn}_r2.fixed.fastq.gz

alignment_dir=${root_dir}/ATAC_seq/alignment_bowtie/
mkdir -p ${alignment_dir}
alignment_fn=${alignment_dir}/${fn}.bowtie2.grch38

samFile=${alignment_fn}.sam
bowtie2 \
  -t \
  -x ${genome_index} \
  -1 ${fastq_dir}/${fn1} \
  -2 ${fastq_dir}/${fn2} \
  -S ${samFile} \
  -p 24

bamFile=${alignment_fn}.sortedByCoord.out.bam

# Keep only reads on autosomal genome
cd ${alignment_dir}
bamFile=${fn}.bowtie2.grch38.sortedByCoord.out.bam

samtools sort ${samFile} > ${bamFile}_intermediate
samtools index ${bamFile}_intermediate
samtools view -@ 24 -b ${bamFile}_intermediate chr{1..22} > ${bamFile}

rm ${samFile}
rm ${bamFile}_intermediate*


