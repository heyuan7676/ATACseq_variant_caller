#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=12:00:00
#SBATCH -p shared

file_dir=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR/ATAC_seq/ftp.sra.ebi.ac.uk/vol1/ERZ683
sample="$1"
sample_dir=${file_dir}/${sample}
#sample=ERZ683839

echo "[INFO]: "$sample

# 0). download the dataset
echo "[INFO]: Downloading..."
bash download_ENA.sh ${sample}

ml samtools

grch37=/work-zfs/abattle4/heyuan/database/GRCh37_reference/GRCh37.p13.genome.fa
cram_fn=`ls ${file_dir}/${sample}/*cram`
donor_ID=`ls ${file_dir}/${sample}/*cram | xargs -n1 basename`
donor_ID=${donor_ID/.cram/}

echo "[INFO]: donor ID - "$donor_ID

# Stop if not in the second batch
N=`cat samples.txt | grep ${donor_ID} | wc -l`
if [ $N -eq 0 ]
then
	echo "[INFO]: donor sample not in the second batch"
	exit
fi

# 1). .cram --> .bam
echo "[INFO]: Converting cram to bam .."
bam_fn=${sample_dir}/${donor_ID}.bam
samtools view -@ 6 -b \
  -T ${grch37} \
  -o ${bam_fn} \
  ${cram_fn}
rm ${cram_fn}

# 2). remove reads mapped to mitochondrial genome (mtDNA)
echo "[INFO]: Removing reads that don't map to autosomal chromosomes..."
bash check_mtDNA.sh ${bam_fn}
samtools view -@ 6 -b ${bam_fn} 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 > ${sample_dir}/${donor_ID}.autosomal.bam
rm ${bam_fn}
bam_fn=${sample_dir}/${donor_ID}.autosomal.bam

# 3). .bam --> paired-ended .fastq files
echo "[INFO]: Converting bam to paired-ended fastq..."
fq_fn=${sample_dir}/${donor_ID}
samtools bam2fq -@ 6 ${bam_fn} > ${fq_fn}.fastq
rm ${bam_fn}

cat ${fq_fn}.fastq | grep '^@.*/1$' -A 3 --no-group-separator > ${fq_fn}_r1.fastq 
cat ${fq_fn}.fastq | grep '^@.*/2$' -A 3 --no-group-separator > ${fq_fn}_r2.fastq
rm ${fq_fn}.fastq

gzip ${fq_fn}_r1.fastq
gzip ${fq_fn}_r2.fastq

# 4). repair the discrepancy in the two paired-end files
# need to provide enough memory, otherwise would crash in the middle and hang in there forever
echo "[INFO]: Repair the discrepancy in the two paired-end files..."
ml bbmap
repair.sh in1=${fq_fn}_r1.fastq.gz in2=${fq_fn}_r2.fastq.gz out1=${fq_fn}_r1.fixed.fastq.gz out2=${fq_fn}_r2.fixed.fastq.gz outsingle=${fq_fn}_singletons.fq.gz
rm ${fq_fn}_r1.fastq.gz
rm ${fq_fn}_r2.fastq.gz

mv ${sample_dir}/${donor_ID}*gz ${file_dir}/