#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH -p shared

ROOT_DIR="/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie"
peak_script_dir=/work-zfs/abattle4/heyuan/Variant_calling/scripts/Peaks

## Call Peaks using MACS2
sbatch ${peak_script_dir}/run_peak.sh -s ${peak_script_dir}/MACS2_call_merge/peak_calling_MACS2.snakefile

bash ${peak_script_dir}/MACS2_call_merge/peak_union_MACS2.sh ${ROOT_DIR}

sbatch ${peak_script_dir}/run_peak.sh -s ${peak_script_dir}/MACS2_call_merge/obtain_peak_matrix_MACS2.snakefile


## Peak correction

PEAK_DIR='${ROOT_DIR}/Peaks_Genrich/combined/'
BAM_DIR='${ROOT_DIR}/first_pass_bqsr'

bash ${peak_script_dir}/run_correct_peaks.sh ${PEAK_DIR} ${BAM_DIR} 
