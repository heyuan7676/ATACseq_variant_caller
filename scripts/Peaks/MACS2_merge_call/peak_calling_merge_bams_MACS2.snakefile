import os
configfile: 'config.yaml'

SAMTOOLS = config['SAMTOOLS']
MACS2 = config['MACS2']

BOWTIE_DIR = config['BOWTIE_DIR']
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks/')
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')

INDIVS = glob_wildcards(os.path.join(DIR_FIRST_PASS, '{indiv}-clean.bam'))
INDIVS = INDIVS[0]
SUFFIX = 'combined_peaks_all.bam'

''' Snakemake rules '''
rule all:
    input:
        #os.path.join(DIR_FIRST_PASS, SUFFIX)
        os.path.join(PEAK_DIR, 'combined_peaks.narrowPeak')

rule merge_bam:
    input:
        expand(os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam'), indiv = INDIVS)
    output:
        os.path.join(DIR_FIRST_PASS, SUFFIX)
    shell:
        '{SAMTOOLS} merge {output} {input}'

rule call_peaks:
    input:
        os.path.join(DIR_FIRST_PASS, SUFFIX)
    output:
        os.path.join(PEAK_DIR, 'combined_peaks.narrowPeak')
    params:
        outdir =  PEAK_DIR,
        prefix =  'combined,
    conda:
        "envs/env_py37.yml"
    shell:
        """
        {MACS2} callpeak -t {input} -f BAM -g hs --nomodel --extsize 147 --max-gap 100 --outdir {params.outdir} -n {params.prefix} -q 0.05
        """


