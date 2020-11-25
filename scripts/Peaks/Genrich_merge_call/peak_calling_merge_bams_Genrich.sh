#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=72:00:00
#SBATCH -p lrgmem

ROOT_DIR="$1"
PEAK_DIR='${ROOT_DIR}/Peaks_Genrich'

DIR_FIRST_PASS=${ROOT_DIR}/first_pass_bqsr
BLACKLIST=${ROOT_DIR}/Peaks/hg38.blacklist.bed

str=${DIR_FIRST_PASS}/`head -n1 test_all.txt`"-clean.sorted.bam"
for s in `cat test_all.txt | sed 1d`
do
	str=${DIR_FIRST_PASS}/${s}"-clean.sorted.bam,"${str}
done

echo $str

signal=${ROOT_DIR}/Peaks_Genrich/combined_peaks.Genrich.narrowPeak
bedgraph=${ROOT_DIR}/Peaks_Genrich/combined_peaks.Genrich.bedgraph

timestamp() {
  echo "Time stamp:" 
  date +"%T" # current time
}
timestamp


/work-zfs/abattle4/heyuan/tools/Genrich-master/Genrich -t ${str} -j -r -o ${signal} -f ${bedgraph} -v -E ${BLACKLIST} -q 0.05


timestamp

