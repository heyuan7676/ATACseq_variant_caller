#!/bin/bash

sample="$1"

# Download metadata

output=records/sample_metadata.tsv

sample=ERR1736146
wget "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${sample}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,base_count,study_title,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true" -O - | head -n1 >  ${output}

for sample in `awk '{print $1}' records/ena_run.tsv | sed 's/"//g' | sed "1d"`
do
	echo $sample
	wget "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${sample}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,base_count,study_title,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true" -O - | sed "1d" >>  ${output}

done


cd /work-zfs/abattle4/heyuan/Variant_calling/datasets/ENA/fastq/
kingfisher get -r ${sample} -m ena-ascp aws-http prefetch
