import os

configfile: 'config.yaml'

'''Load from config file'''

BOWTIE_DIR = config['BOWTIE_DIR']
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')
PEAK_DIR_MACS2 = os.path.join(BOWTIE_DIR, 'Peaks_MACS2')
PEAK_DIR_Genrich = os.path.join(BOWTIE_DIR, 'Peaks_Genrich')

os.makedirs(PEAK_DIR_MACS2, exist_ok = True)
os.makedirs(PEAK_DIR_Genrich, exist_ok = True)

INDIVS = glob_wildcards(os.path.join(DIR_FIRST_PASS, '{indiv}-clean.bam'))
INDIVS = INDIVS[0]
INDIVS = [ id for id in INDIVS if '/' not in id ]
INDIVS.sort()

BEDTOOLS = config['BEDTOOLS']
SAMTOOLS = config['SAMTOOLS']
MACS2 = config['MACS2']
GENRICH = config['GENRICH']

## These need to after the global variables in order to use them 
include : '../scripts/Peaks/peak_calling.snakefile'


''' Snakemake rules '''
rule all:
    input:
        expand(os.path.join(PEAK_DIR_MACS2, '{indiv}' + '_peaks.narrowPeak'), indiv = INDIVS),
        #expand(os.path.join(PEAK_DIR_Genrich, '{indiv}' + '_peaks.narrowPeak'), indiv = INDIVS)
