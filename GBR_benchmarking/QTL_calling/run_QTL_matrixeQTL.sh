#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=6:00:00
#SBATCH -p lrgmem

chromosome="$1"
minDP=3
ml python/3.7


#python QTL_calling_prepare_for_matrixeQTL.py ${minDP} ${chromosome}
python temp.py ${minDP} ${chromosome}

ml gcc/5.5.0
ml R/3.5.1

#for dist in 0 1000 10000 100000 1000000
for dist in 0 1000 10000 100000 1000000
do
	#for method in WGS GC Imputation Integration
        for method in Integration
	do
		Rscript matrix_eQTL.R ${chromosome} ${method} ${dist}
	done
done


