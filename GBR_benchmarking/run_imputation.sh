#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=shared

#INDIV="$1" snakemake -s Snakefile --cores 24

source ~/.bash_profile

ml gcc/5.5.0
ml picard
ml bcftools
ml vcftools
ml htslib

source activate snakemake

## Run one snakemake job
#INDIV="$1" snakemake -s Snakefile --cores 24


## Run multiple jobs in parallel
snakemake --use-conda --jobs 20 \
    --cluster "sbatch --ntasks=12 --time=6:00:00 --partition shared"\
    --rerun-incomplete \
    --keep-going \
    --latency-wait 5 \
    --printshellcmds \
    $@

# Note: need to specify conda environment in rule

