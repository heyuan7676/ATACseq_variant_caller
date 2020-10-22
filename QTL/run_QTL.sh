#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=24:00:00
#SBATCH --partition shared

ml python/2.7
python QTL_calling.py "$1" "$2" "$3"
