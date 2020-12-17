#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH -p shared


# Imputation - step 1

# Imputation - step 2
sbatch run_imputation.sh -s /work-zfs/abattle4/heyuan/Variant_calling/scripts/Genotype_imputation/obtain_genotype_imputed.snakefile
