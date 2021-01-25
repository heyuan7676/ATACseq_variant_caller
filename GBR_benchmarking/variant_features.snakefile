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

PEAK_DIR_MACS2 = os.path.join(BOWTIE_DIR, 'Peaks_MACS2')
PEAK_DIR_Genrich = os.path.join(BOWTIE_DIR, 'Peaks_Genrich')

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
minDP_arr = ['2','3', '4', '6']

CHROM = config['CHROM']
wildcard_constraints:
    chr="\d+"



## These need to after the global variables in order to use them 
include : '../scripts/Genotype_calling/collect_variant_features.snakefile'

output_suffix = 'recode'
''' Snakemake rules '''
rule all:
    input:
        expand(os.path.join(VCF_DIR, 'minDP{minDP}', "{indiv}.filtered.recode.INFO.txt"),  minDP = minDP_arr, indiv = INDIVS),
        expand(os.path.join(VCF_DIR, "minDP{minDP}", '{indiv}' + ".filtered.minDP{minDP}.recode.variants.toMACS2Peaks.txt"), minDP = minDP_arr, indiv = INDIVS),
        expand(os.path.join(IMPUTE_DIR, "minDP{minDP}", '{indiv}' + ".filtered.minDP{minDP}.imputed.variants.toMACS2Peaks.txt"), minDP = minDP_arr, indiv = INDIVS),
