#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=6:00:00
#SBATCH -p shared

chromosome="$1"
ml python/3.7


python QTL_calling_prepare_for_matrixeQTL.py 3 ${chromosome}

ml gcc/5.5.0
ml R/3.5.1

for dist in 0 1000 10000 100000 1000000
do
	for method in GC Imputation Integration
	do
		Rscript matrix_eQTL.R ${chromosome} ${method} ${dist}
	done
done
