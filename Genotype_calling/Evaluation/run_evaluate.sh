#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH -p skylake

sample="$1"
#fn=performance/combined/called_in_both_${sample}_allSNPs.txt
#if [ ! -f "$fn" ]
#then
	ml python/2.7
	python Evaluation_metrics_combined.py ${sample}
#fi
