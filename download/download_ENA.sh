#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH -p skylake

sample="$1"
echo $sample
f1=`cat SRA_file_reports.txt | grep $sample | awk -F'	' '{print $7}' | cut -d';' -f1`
f2=`cat SRA_file_reports.txt | grep $sample | awk -F'	' '{print $7}' | cut -d';' -f2`


cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/ftp.sra.ebi.ac.uk/vol1/ERZ683
mkdir -p ${sample}
cd ${sample}

wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ683/${sample}/${f1}
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ683/${sample}/${f2}

