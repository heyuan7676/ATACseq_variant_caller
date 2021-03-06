#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH -p skylake

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
snakemake --use-conda --jobs 40 \
    --cluster "sbatch --ntasks=4 --time=24:00:00 --partition skylake" \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 5 \
    --printshellcmds \
    $@

# Note: need to specify conda environment in rule


## failed
#~/miniconda3/envs/snakemake/bin/snakemake -s Snakefile \
#          --profile /home-4/yhe23@jhu.edu/.config/snakemake/slurm \
#          --rerun-incomplete \
#          --use-conda
# Problem: the submitted job can't find the snakemake module
# Error: No module named snakemake.utils


