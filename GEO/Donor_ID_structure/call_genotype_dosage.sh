#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH -p skylake

ml python/2.7
root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Subsampling50
for minDP in 3 4 5 6 8 10
do
	for s in `cat Test_samples.txt`
	do
		fn=${root_dir}/VCF_files/minDP${minDP}/${s}.filtered.minDP${minDP}.recode.dosage_genotype.bed
		echo "minDP${minDP}: sample $s"
		python ../../scripts/Genotype_calling/Derive_dosage_mll.py ${minDP} ${s} ${root_dir}
	done
done
