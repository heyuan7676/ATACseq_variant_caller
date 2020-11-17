#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --time=6:00:00
#SBATCH -p shared

ml python/2.7

chromosome="$1"
python QTL_calling_run.py ${chromosome} 0 True
python QTL_calling_run.py ${chromosome} 0 False
