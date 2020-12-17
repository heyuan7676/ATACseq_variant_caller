#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=6:00:00

#INDIV="$1" snakemake -s Snakefile --cores 24


ml bedtools
ml samtools
source activate snakemake

## Run one snakemake job
#INDIV="$1" snakemake -s Snakefile --cores 24


## Run multiple jobs in parallel
snakemake --use-conda --jobs 2 \
    --cluster "sbatch --ntasks=24 --time=24:00:00 --partition lrgmem"\
    --rerun-incomplete \
    --keep-going \
    --latency-wait 5 \
    --printshellcmds \
    $@

# Note: need to specify conda environment in rule

