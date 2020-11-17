import os

configfile: 'config.yaml'

'''Load from config file'''

PICARD = config['PICARD']
MACS2 = config['MACS2']
SAMTOOLS = config['SAMTOOLS']

BOWTIE_DIR = config['BOWTIE_DIR']
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks')

INDIVS = glob_wildcards(os.path.join(DIR_FIRST_PASS, '{indiv}-clean.bam'))
INDIVS = INDIVS[0]

''' Snakemake rules '''
rule all:
    input:
        expand(os.path.join(PEAK_DIR, '{indiv}' + '_peaks.narrowPeak'), indiv = INDIVS)


'''Call peaks using MACS2'''
rule peak_calling:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam')
    output:
        os.path.join(PEAK_DIR, '{indiv}' + '_peaks.narrowPeak')
    params:
        outdir =  PEAK_DIR,
        prefix =  '{indiv}',
    conda:
        "envs/env_py37.yml"
    shell:
        """
        {MACS2} callpeak -t {input} -f BAM -g hs --shift -75 --extsize 150 --max-gap 100 --outdir {params.outdir} -n {params.prefix} -q 0.05
        """


