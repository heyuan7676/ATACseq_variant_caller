#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH -p shared

source activate snakemake
#INDIV="$1" snakemake -s Snakefile --cores 24

ml picard
ml samtools
ml gatk
ml bowtie2
ml bcftools
ml vcftools

## Run one snakemake job
#INDIV="$1" snakemake -s Snakefile --cores 24


## Run multiple jobs in parallel
snakemake --use-conda --jobs 20 \
    --cluster "sbatch --ntasks=1 --time=6:00:00 --partition shared" \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
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


