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
INDIVS.sort()
INDIVS = INDIVS[3:5]

''' Snakemake rules '''
rule all:
    input:
        expand(os.path.join(PEAK_DIR, '{indiv}' + '_useBED_peaks.narrowPeak'), indiv = INDIVS)



rule bam_to_bed:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam')
    output:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bed')
    shell:
        """
        bedtools bamtobed -i {input} > {output}
        """


'''Call peaks using MACS2'''
rule peak_calling:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bed')
    output:
        os.path.join(PEAK_DIR, '{indiv}' + '_useBED_peaks.narrowPeak')
    params:
        outdir =  PEAK_DIR,
        prefix =  '{indiv}_useBED',
    conda:
        "envs/env_py37.yml"
    shell:
        """
        {MACS2} callpeak -t {input} -f BED -g hs --nomodel --shift -100 --extsize 200 --max-gap 100 --outdir {params.outdir} -n {params.prefix} -q 0.05
        """


