#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=12:00:00
#SBATCH -p parallel

file_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/ftp.sra.ebi.ac.uk/vol1/ERZ683
sample="$1"
sample_dir=${file_dir}/${sample}
#sample=ERZ683839

f1=`cat SRA_file_reports.txt | grep $sample | awk -F'	' '{print $7}' | cut -d';' -f1 | xargs -n1 basename`
f2=`cat SRA_file_reports.txt | grep $sample | awk -F'	' '{print $7}' | cut -d';' -f2 | xargs -n1 basename`

if [[ $f1 == *cram ]]
then
	donor_ID=${f1/.cram/}
else
	donor_ID=${f2/.cram/}
fi

echo "[INFO]: Sample - "$sample" ; "donor ID - $donor_ID

outputFN=${file_dir}/${donor_ID}_r1.fixed.fastq.gz


# 0). download the dataset
#echo "[INFO]: Downloading..."
bash download_ENA.sh ${sample}

grch37=/work-zfs/abattle4/heyuan/database/GRCh37_reference/GRCh37.p13.genome.fa
# if download fails, exist
N=`ls ${file_dir}/${sample}/* | wc -l`
if [ $N -lt 2 ]
then
       echo "[INFO]: download fails, exit"
       exit
fi

# Stop if not in the second batch
#N=`cat samples.txt | grep ${donor_ID} | wc -l`
#if [ $N -eq 0 ]
#then
#	echo "[INFO]: donor sample not in the second batch"
#	exit
#fi


ml samtools
# 1). .cram --> .bam
echo "[INFO]: Converting cram to bam .."
bam_fn=${sample_dir}/${donor_ID}.bam
cram_fn=${sample_dir}/${donor_ID}.cram
samtools view -@ 24 -b -T ${grch37} -o ${bam_fn} ${cram_fn}
#rm ${cram_fn}

# 2). remove reads mapped to mitochondrial genome (mtDNA)
echo "[INFO]: Removing reads that don't map to autosomal chromosomes..."
bash check_mtDNA.sh ${bam_fn}
samtools view -@ 24 -b ${bam_fn} 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 > ${sample_dir}/${donor_ID}.autosomal.bam
#rm ${bam_fn}

# 3).bam --> fastq files
echo "[INFO]: Converting bam to fastq..."
bam_auto_fn=${sample_dir}/${donor_ID}.autosomal.bam
fq_fn=${sample_dir}/${donor_ID}
samtools bam2fq -@ 24 ${bam_auto_fn} > ${fq_fn}.fastq
#rm ${bam_auto_fn}

# 4).remove reads with duplicated names in the fastq file
cat ${fq_fn}.fastq | seqkit rmdup -n -o ${fq_fn}.clean.fastq.gz -d ${fq_fn}.duplicated.fastq.gz -D ${fq_fn}.duplicated.txt
echo "[INFO]: ${fq_fn}.fastq"
# rm ${fq_fn}.fastq

# 5).split in the paired-ended files
zcat ${fq_fn}.clean.fastq.gz | grep '^@.*/1$' -A 3 --no-group-separator > ${fq_fn}_r1.fastq 
zcat ${fq_fn}.clean.fastq.gz | grep '^@.*/2$' -A 3 --no-group-separator > ${fq_fn}_r2.fastq


# 6). for samples in second batch
gzip ${fq_fn}_r1.fastq
gzip ${fq_fn}_r2.fastq
echo "[INFO]: ${fq_fn}.clean.fastq.gzq"
# rm ${fq_fn}.clean.fastq.gz


# 6). change quanlity encoding for the samples in first batch
#reformat.sh in=${fq_fn}_r1.fastq out=${fq_fn}_r1.qin.fastq.gz qin=33 qout=64 tossbrokenreads
#reformat.sh in=${fq_fn}_r2.fastq out=${fq_fn}_r2.qin.fastq.gz qin=33 qout=64 tossbrokenreads
#mv ${fq_fn}_r1.qin.fastq.gz ${fq_fn}_r1.fastq.gz
#mv ${fq_fn}_r2.qin.fastq.gz ${fq_fn}_r2.fastq.gz


# 6). repair the discrepancy in the two paired-end files
# need to provide enough memory, otherwise would crash in the middle and hang in there forever
# use the most updated bbmap version!

# if the sample was done in the first batch, needs to change the encoding
#reformat.sh in=${fq_fn}_r1.fastq.gz out=${fq_fn}_r1.qin.fastq.gz qin=33 qout=64
#reformat.sh in=${fq_fn}_r2.fastq.gz out=${fq_fn}_r2.qin.fastq.gz qin=33 qout=64

echo "[INFO]: Repair the discrepancy in the two paired-end files..."
repair.sh in1=${fq_fn}_r1.fastq.gz in2=${fq_fn}_r2.fastq.gz out1=${fq_fn}_r1.fixed.fastq.gz out2=${fq_fn}_r2.fixed.fastq.gz outsingle=${fq_fn}_singletons.fq.gz
echo "[INFO]: ${fq_fn}_r1.fastq.gz"
# rm ${fq_fn}_r1.fastq.gz
# rm ${fq_fn}_r2.fastq.gz

#mv ${fq_fn}_r1.fixed.fastq.gz ${file_dir}
#mv ${fq_fn}_r2.fixed.fastq.gz ${file_dir}
