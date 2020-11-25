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
ml bedtools

## Run one snakemake job
#INDIV="$1" snakemake -s Snakefile --cores 24


## Run multiple jobs in parallel
snakemake --use-conda --jobs 40 \
    --cluster "sbatch --ntasks=2 --time=24:00:00 --partition shared" \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 5 \
    --printshellcmds \
    $@

