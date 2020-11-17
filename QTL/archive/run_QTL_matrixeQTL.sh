#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=6:00:00
#SBATCH -p shared

chromosome="$1"

#ml R/3.5.1

#for dis in 0 1000
#do
#	Rscript matrix_eQTL.R ${chromosome} "" ${dis}
#	Rscript matrix_eQTL.R ${chromosome} "_realGT" ${dis}
#	Rscript matrix_eQTL.R ${chromosome} "_realGT_all" ${dis}
#done
