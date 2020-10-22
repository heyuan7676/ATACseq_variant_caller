#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition shared

ml python/2.7
python Evaluation_metrics_QTLs.py "$1" "$2"
