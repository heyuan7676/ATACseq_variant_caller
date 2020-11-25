#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH -p shared



ml python/2.7
peak_dir="$1"
bam_dir="$2"

#PEAK_DIR='/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks_Genrich/combined/'
#BAM_DIR='/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/first_pass_bqsr'


script_dir=/work-zfs/abattle4/heyuan/Variant_calling/scripts/Peaks
python ${script_dir}/Correct_peaks.py ${peak_dir} ${bam_dir}
