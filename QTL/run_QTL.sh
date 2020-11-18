#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=6:00:00
#SBATCH -p skylake

ml python/2.7

chromosome="$1"

for peak_calling in macs2 Genrich
do
	python QTL_calling_run.py ${chromosome} 0 ${peak_calling}
done

