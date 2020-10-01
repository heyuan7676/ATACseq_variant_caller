#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=24:00:00
#SBATCH -p shared

ml bowtie2
ml samtools

genome_index=/work-zfs/abattle4/heyuan/database/bowtie2/grch37/GRCh37

fn=HG00247
echo $fn

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/
fastq_dir=${root_dir}/ATAC_seq/ftp.sra.ebi.ac.uk/vol1/ERZ683/
fn1=${fn}_r1.fastq.gz
fn2=${fn}_r2.fastq.gz

alignment_dir=${root_dir}/ATAC_seq/alignment_bowtie/
mkdir -p ${alignment_dir}
alignment_fn=${alignment_dir}/${fn}.bowtie2.grch37

samFile=${alignment_fn}.sam
bowtie2 \
  -x ${genome_index} \
  -1 ${fastq_dir}/${fn1} \
  -2 ${fastq_dir}/${fn2} \
  -S ${samFile} \
  -p 24


bamFile=${alignment_fn}.sortedByCoord.out.bam
samtools view -bS ${samFile} | samtools sort -o ${bamFile}


