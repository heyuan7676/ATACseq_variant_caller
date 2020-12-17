#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH -p shared

#INDIV="$1" snakemake -s Snakefile --cores 24


ml bedtools

source activate snakemake

## Run one snakemake job
#INDIV="$1" snakemake -s Snakefile --cores 24


## Run multiple jobs in parallel
snakemake --use-conda --jobs 5 \
    --cluster "sbatch --ntasks=4 --time=6:00:00 --partition shared"\
    --rerun-incomplete \
    --keep-going \
    --latency-wait 5 \
    --printshellcmds \
    $@

# Note: need to specify conda environment in rule

