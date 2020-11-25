#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH -p shared

file_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/ftp.sra.ebi.ac.uk/vol1/ERZ683
sample=ERZ683899
sample_dir=${file_dir}/${sample}

f1=`cat SRA_file_reports.txt | grep $sample | awk -F'	' '{print $7}' | cut -d';' -f1 | xargs -n1 basename`
f2=`cat SRA_file_reports.txt | grep $sample | awk -F'	' '{print $7}' | cut -d';' -f2 | xargs -n1 basename`
if [[ $f1 == *cram ]]
then
        donor_ID=${f1/.cram/}
else
        donor_ID=${f2/.cram/}
fi

echo $f1
echo $f2


# 0). download the dataset
#echo "[INFO]: Downloading..."
#bash download_ENA.sh ${sample}

grch37=/work-zfs/abattle4/heyuan/database/GRCh37_reference/GRCh37.p13.genome.fa
cram_fn=${sample_dir}/${donor_ID}.cram

ml samtools
samtools view -b  -T ${grch37} -o ${cram_fn%.cram}.bam ${cram_fn}

