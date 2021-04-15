#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH -p skylake

for peak_calling in MACS2 MACS2/combined Genrich Genrich/combined
do

for method in VCF_files Imputation Integration
do
outdir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/fastQTL/${method}/${peak_calling}

for cisDis in 100 300 500 1000 10000 100000 1000000
do
    outfn=${outdir}/chr*.window${cisDis}.fastq.permutation.results
    N=`ls ${outdir}/chr*.window${cisDis}.fastq.permutation.results | wc -l`
    if [[ "$N" != 22 ]]
    then
        echo "FastQTL not finished: ", ${method}, ${peak_calling}, ${cisDis}
        continue
    fi
#    FILE=${outdir}/window${cisDis}.fastq.permutation.results.BH.txt
#    if test -f "$FILE"; then
#        echo "$FILE exists."
#	continue
#    fi

    echo $method, $peak_calling, $cisDis
    cat ${outfn} > ${outdir}/window${cisDis}.fastq.permutation.results.txt
    Rscript collect_fastQTL.R ${peak_calling} ${method} ${cisDis}

done


done

done
