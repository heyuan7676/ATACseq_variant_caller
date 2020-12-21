import os

configfile: 'config.yaml'

'''Load from config file'''

PICARD = config['PICARD']
GATK = config['GATK']
SAMTOOLS = config['SAMTOOLS']
BCFTOOLS = config['BCFTOOLS']
VCFTOOLS = config['VCFTOOLS']
BEDTOOLS = config['BEDTOOLS']
MACS2 = config['MACS2']

DIR = config['DIR']

INDIVS = []
fn = open('samples.txt', 'r')
for line in fn.readlines():
    INDIVS.append(line.rstrip())
INDIVS = INDIVS[:10]


FQ_DIR = config['FQ_DIR']
BOWTIE_DIR = config['BOWTIE_DIR']
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks')
VCF_DIR = os.path.join(BOWTIE_DIR, 'VCF_files') 
GENOTYPE_DIR = os.path.join(BOWTIE_DIR, 'Called_GT')
TMP_DIR = config['TMP_DIR']

os.makedirs(FQ_DIR, exist_ok = True)
os.makedirs(BOWTIE_DIR, exist_ok = True)
os.makedirs(DIR_FIRST_PASS, exist_ok = True)
os.makedirs(PEAK_DIR, exist_ok = True)
os.makedirs(VCF_DIR, exist_ok = True)
os.makedirs(GENOTYPE_DIR, exist_ok = True)
os.makedirs(TMP_DIR, exist_ok = True)

# Bowtie 2 genome
BOWTIE_GENOME_INDEX = config['BOWTIE_GENOME_INDEX']
SUFFIX = '.bowtie2.grch38.sortedByCoord.out'

GENOME = config['GENOME']
GENOME_STAR = config['GENOME_STAR']
VCFFN = config['VCFFN']

THREADS = config['THREADS']
minDP_arr = config['minDP_arr']

for dp in minDP_arr:
    os.makedirs(os.path.join(VCF_DIR, 'minDP' + dp), exist_ok = True)
    os.makedirs(os.path.join(GENOTYPE_DIR, 'minDP'+ dp), exist_ok = True)
    os.makedirs(os.path.join(VCF_DIR, 'minDP' + dp, 'GRCh37'), exist_ok = True)
    for s in INDIVS:
        os.makedirs(os.path.join(VCF_DIR, 'minDP' + dp, 'GRCh37', s), exist_ok = True)
        os.makedirs(os.path.join(BOWTIE_DIR, 'Imputation', 'minDP' + dp, s), exist_ok = True)


CHROM = config['CHROM']
wildcard_constraints:
    chr="\d+"



## These need to after the global variables in order to use them 
include : '../scripts/Genotype_calling/alignment.snakefile'
include : '../scripts/Genotype_calling/processing_QC.snakefile'
include : '../scripts/Genotype_calling/variant_calling.snakefile'
include : '../scripts/Genotype_calling/obtain_genotype.snakefile'
include : '../scripts/Genotype_calling/genotype_imputation.snakefile'
#include : '../scripts/Genotype_calling/obtain_genotype_imputed.snakefile'


''' Snakemake rules '''
rule all:
    input:
        expand(os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.filtered.minDP' + '{minDP}' + '.recode.vcf.gz'), indiv = INDIVS, minDP = minDP_arr),
        expand(os.path.join(GENOTYPE_DIR, 'minDP{minDP}', '{indiv}.filtered.genotype.minDP' + '{minDP}' + '.txt'), indiv = INDIVS, minDP = minDP_arr),
        expand(os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.formatted.vcf"), indiv = INDIVS),
        expand(os.path.join(BOWTIE_DIR, 'Imputation', 'minDP' + "{minDP}", "{indiv}", "{indiv}.imputed.dosage.GRCh38.bed"), minDP = minDP_arr, indiv = INDIVS)
        #expand(os.path.join(GENOTYPE_DIR, "minDP{minDP}", 'union-SNPs_minDP{minDP}.bed'), minDP = minDP_arr),
        #expand(os.path.join(GENOTYPE_DIR, "minDP{minDP}", "gt_by_sample_matrix_chr{chr}.txt"), chr = CHROM, minDP = minDP_arr),
        #expand(os.path.join(VCF_DIR, "gt_info_by_sample_matrix_chr{chr}_atac.txt"), chr = CHROM)
        #os.path.join(IMPUTATION_DIR,  "union-SNPs_imputed.bed"),
        #expand(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.filtered.genotype.txt_matrix"), indiv = INDIVS),
        #expand(os.path.join(IMPUTATION_DIR, "gt_by_sample_matrix_chr{chr}.txt"), chr = CHROM)

