#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00

source activate snakemake
#INDIV="$1" snakemake -s Snakefile --cores 24

ml samtools
ml bedtools

## Run one snakemake job
INDIV="$1" snakemake --use-conda --printshellcmds -s Snakefile --cores 1


## Run multiple jobs in parallel
#snakemake --use-conda --jobs 20 \
#    --cluster "sbatch --ntasks=1 --time=6:00:00" \
#    --rerun-incomplete \
#    --keep-going \
#    --latency-wait 10 \
#    --printshellcmds \
#    $@

# Note: need to specify conda environment in rule


