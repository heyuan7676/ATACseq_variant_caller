#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=20:00:00
#SBATCH -p skylake

# Access GEO
# Obtained the metadata and accession list from SRA Run Selector
# HSC_AccessList.txt  HSC_metadata.txt

# HSC: GSE74912, GSE74246
# ImmGen: GSE100738, GSE109125

AccessList_Fn=Yoruba/AccessList.txt


# Download the SRA files 
Dataset=Yoruba
DataType=DNase_seq
save_folder=/work-zfs/abattle4/heyuan/Variant_calling/benchmarking/datasets/${Dataset}/${DataType}/download

mkdir -p ${save_folder}

ml sra-tools
prefetch --option-file ${AccessList_Fn} -O ${save_folder}


cd ${save_folder}
for fn in `ls *sra`
do
	echo $fn
	fastq-dump --split-files ${fn}
done
