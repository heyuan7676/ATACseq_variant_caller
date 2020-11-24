#B/bin/bash

sample="$1"

bash split_chrs.sh ${sample}
bash liftHg38ToHg19.sh ${sample}
bash imputation.sh ${sample}
bash obtain_genotype.sh ${sample}
