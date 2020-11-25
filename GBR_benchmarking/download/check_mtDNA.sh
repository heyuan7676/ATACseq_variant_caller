#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --time=2:00:00
#SBATCH -p debug

#!/bin/bash

#### Calculate percentage of reads mapped to mitochondrial genome (mtDNA) using SAMtools idxstats
#### Can be useful for ATAC-seq data. Requires an indexed BAM file:

## Check if index is present. If not, create it:
if [[ ! -e ${1}.bai ]];
  then
  echo '[INFO]: File does not seem to be indexed. Indexing now:'
  samtools index $1
fi

## Calculate %mtDNA:
mtReads=$(samtools idxstats $1 | grep 'MT' | cut -f 3)
totalReads=$(samtools idxstats $1 | awk '{SUM += $3} END {print SUM}')

echo '[INFO]: mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'
