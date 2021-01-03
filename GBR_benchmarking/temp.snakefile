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

FQ_DIR = config['FQ_DIR']
BOWTIE_DIR = config['BOWTIE_DIR']
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks')
VCF_DIR = os.path.join(BOWTIE_DIR, 'VCF_files') 
IMPUTE_DIR = os.path.join(BOWTIE_DIR, 'Imputation')
INTERGRATED_DIR = os.path.join(BOWTIE_DIR, 'Integration')
GENOTYPE_DIR = os.path.join(BOWTIE_DIR, 'Called_GT')
TMP_DIR = config['TMP_DIR']


os.makedirs(FQ_DIR, exist_ok = True)
os.makedirs(BOWTIE_DIR, exist_ok = True)
os.makedirs(DIR_FIRST_PASS, exist_ok = True)
os.makedirs(PEAK_DIR, exist_ok = True)
os.makedirs(VCF_DIR, exist_ok = True)
os.makedirs(IMPUTE_DIR, exist_ok = True)
os.makedirs(INTERGRATED_DIR, exist_ok = True)
os.makedirs(GENOTYPE_DIR, exist_ok = True)
os.makedirs(TMP_DIR, exist_ok = True)

# Bowtie 2 genome
BOWTIE_GENOME_INDEX = config['BOWTIE_GENOME_INDEX']
SUFFIX = '.bowtie2.grch38.sortedByCoord.out'
GENOME = config['GENOME']
VCFFN = config['VCFFN']

THREADS = config['THREADS']
minDP_arr = config['minDP_arr']
minDP_arr = ['4']
for dp in minDP_arr:
    os.makedirs(os.path.join(VCF_DIR, 'minDP' + dp), exist_ok = True)
    os.makedirs(os.path.join(GENOTYPE_DIR, 'minDP'+ dp), exist_ok = True)
    os.makedirs(os.path.join(VCF_DIR, 'minDP' + dp, 'GRCh37'), exist_ok = True)
    for s in INDIVS:
        os.makedirs(os.path.join(VCF_DIR, 'minDP' + dp, 'GRCh37', s), exist_ok = True)
        os.makedirs(os.path.join(IMPUTE_DIR, 'minDP' + dp, s), exist_ok = True)


CHROM = config['CHROM']
wildcard_constraints:
    chr="\d+"



## These need to after the global variables in order to use them 
include : '../scripts/Genotype_calling/alignment.snakefile'
include : '../scripts/Genotype_calling/processing_QC.snakefile'
include : '../scripts/Genotype_calling/variant_calling.snakefile'
include : '../scripts/Genotype_calling/genotype_imputation.snakefile'
include : '../scripts/Genotype_calling/obtain_genotype_called.snakefile'
include : '../scripts/Genotype_calling/obtain_genotype_imputed.snakefile'
include : '../scripts/Genotype_calling/obtain_genotype_integrated.snakefile'

output_suffix_called = 'recode'
output_suffix_integrated = 'with_Inconsistent'
''' Snakemake rules '''
rule all:
    input:
        expand(os.path.join(IMPUTE_DIR, 'minDP' + "{minDP}", "dosage_by_sample_matrix_chr{chr}.txt"), minDP = minDP_arr, chr = CHROM),
        expand(os.path.join(VCF_DIR, 'minDP' + "{minDP}", "dosage_by_sample_matrix_chr{chr}." + output_suffix_called + ".txt"), minDP = minDP_arr, chr = CHROM),
        expand(os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", "dosage_by_sample_matrix_chr{chr}." + output_suffix_integrated+ ".txt"), minDP = minDP_arr, chr = CHROM)
