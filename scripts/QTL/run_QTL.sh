#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH -p shared

ml python/2.7

chromosome="$1"
window="$2"

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq
for peak_calling in Genrich_combined
do
	python QTL_calling_run.py ${root_dir} ${chromosome} ${window} ${peak_calling}
done

