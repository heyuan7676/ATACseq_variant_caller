#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=10:00:00




wget -O ENCLB544CVI_L1_R1.fa.gz https://www.encodeproject.org/files/ENCFF405CPB/@@download/ENCFF405CPB.fastq.gz
wget -O ENCLB544CVI_L1_R2.fa.gz https://www.encodeproject.org/files/ENCFF846EWZ/@@download/ENCFF846EWZ.fastq.gz
wget -O ENCLB544CVI_L2_R1.fa.gz https://www.encodeproject.org/files/ENCFF521XIS/@@download/ENCFF521XIS.fastq.gz
wget -O ENCLB544CVI_L2_R2.fa.gz https://www.encodeproject.org/files/ENCFF494CZA/@@download/ENCFF494CZA.fastq.gz
wget -O ENCLB544CVI_L3_R1.fa.gz https://www.encodeproject.org/files/ENCFF618XNK/@@download/ENCFF618XNK.fastq.gz
wget -O ENCLB544CVI_L3_R2.fa.gz https://www.encodeproject.org/files/ENCFF673BIU/@@download/ENCFF673BIU.fastq.gz
wget -O ENCLB358EUZ_L1_R1.fa.gz https://www.encodeproject.org/files/ENCFF313TQJ/@@download/ENCFF313TQJ.fastq.gz
wget -O ENCLB358EUZ_L1_R2.fa.gz https://www.encodeproject.org/files/ENCFF413FVI/@@download/ENCFF413FVI.fastq.gz
wget -O ENCLB358EUZ_L2_R1.fa.gz https://www.encodeproject.org/files/ENCFF400MSN/@@download/ENCFF400MSN.fastq.gz
wget -O ENCLB358EUZ_L2_R2.fa.gz https://www.encodeproject.org/files/ENCFF349AHE/@@download/ENCFF349AHE.fastq.gz
