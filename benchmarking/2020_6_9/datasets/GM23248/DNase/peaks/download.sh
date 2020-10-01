#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=10:00:00


wget -O ENCLB544CVI.peaks.bed.gz https://www.encodeproject.org/files/ENCFF767NOX/@@download/ENCFF767NOX.bed.gz
wget -O ENCLB358EUZ.peaks.bed.gz https://www.encodeproject.org/files/ENCFF285EBH/@@download/ENCFF285EBH.bed.gz


gunzip *gz
for x in *bed
do
	sed 's/^chr//g' ${x} > ${x}_temp
	mv ${x}_temp ${x}
done
