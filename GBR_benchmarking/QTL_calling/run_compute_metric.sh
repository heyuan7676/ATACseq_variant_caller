#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH -p lrgmem

ml python/2.7


peak_calling="$1"
window="$2"
python Compute_metrics.py ${peak_calling} ${window}

