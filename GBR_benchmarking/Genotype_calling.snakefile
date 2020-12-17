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
fn = open('test.txt', 'r')
for line in fn.readlines():
    INDIVS.append(line.rstrip())

FQ_DIR = config['FQ_DIR']
BOWTIE_DIR = config['BOWTIE_DIR']
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks')
VCF_DIR = os.path.join(BOWTIE_DIR, 'VCF_files') 
GENOTYPE_DIR = os.path.join(BOWTIE_DIR, 'Called_GT')
TMP_DIR = config['TMP_DIR']

# Bowtie 2 genome
BOWTIE_GENOME_INDEX = config['BOWTIE_GENOME_INDEX']
SUFFIX = '.bowtie2.grch38.sortedByCoord.out'

GENOME = config['GENOME']
GENOME_STAR = config['GENOME_STAR']
VCFFN = config['VCFFN']

THREADS = config['THREADS']
minDP_arr = config['minDP_arr']
GQ_arr = config['GQ_arr']

## These need to after the global variables in order to use them 
include : '../scripts/alignment.snakefile'
include : '../scripts/processing_QC.snakefile'
include : '../scripts/variant_calling.snakefile'



''' Snakemake rules '''
rule all:
    input:
        #expand(os.path.join(BOWTIE_DIR, '{indiv}' + '.bowtie2.grch38.sam'), indiv = INDIVS)
        #expand(os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.sorted.bam.bai'), indiv = INDIVS)
        expand(os.path.join(GENOTYPE_DIR, 'minDP{minDP}', '{indiv}.filtered.genotype.minDP' + '{minDP}' + '.txt'), indiv = INDIVS, minDP = minDP_arr),
        expand(os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.formatted.vcf"), indiv = INDIVS)
        #expand(os.path.join(GENOTYPE_DIR, '{indiv}.filtered.genotype.GQ' + '{GQ}' + '.txt'), indiv = INDIVS, GQ = GQ_arr),
        #expand(os.path.join(VCF_DIR, '{indiv}.filtered.recode.INFO.vcf'), indiv = INDIVS)
