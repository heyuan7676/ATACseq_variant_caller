#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH -p skylake

source ~/.bash_profile

ml gcc/5.5.0
ml picard
ml samtools
ml gatk
ml bowtie2
ml bcftools
ml vcftools
ml bedtools
ml htslib

source activate snakemake

## Run multiple jobs in parallel
snakemake --use-conda --jobs 50 \
    --cluster "sbatch --ntasks=6 --time=6:00:00 --partition skylake" \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 10 \
    --printshellcmds \
    $@


#ml python/2.7
#for minDP in {2..10}
#do
#	for s in `cat samples.txt`
#	do
#		fn=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/VCF_files/minDP${minDP}/${s}.filtered.minDP${minDP}.recode.dosage_genotype.bed
#		echo "minDP${minDP}: sample $s"
#		python ../scripts/Genotype_calling/Derive_dosage_mll.py ${minDP} ${s}
#	done
#done
