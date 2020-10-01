#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH -p skylake


cd /work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/GBR

mkdir -p ATAC_seq
cd ATAC_seq

sample="$1"
echo $sample
wget -r -l3 ftp://ftp.sra.ebi.ac.uk/vol1/ERZ683/${sample}/

