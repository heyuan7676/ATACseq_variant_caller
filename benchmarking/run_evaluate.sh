#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH -p skylake


ml python/2.7
python Evaluation_metrics_Run.py 2
