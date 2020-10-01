#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00


ml fastqc
cd /work-zfs/abattle4/heyuan/ASI/GM23248/DNase/faFiles
for fn in *R1.fa.gz
do
	testFn=${fn%.gz}_fastqc
	if [[ -f "$testFn" ]]
	then
		continue
	fi
	fastqc ${fn}
	unzip ${fn%.gz}_fastqc.zip
	paste <(echo ${fn%.fa.gz}) <(cat ${fn%.gz}_fastqc/fastqc_data.txt | grep "Total Sequence" | awk '{print $3}')
done
