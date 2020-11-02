import os
configfile: 'config.yaml'

BOWTIE_DIR = config['BOWTIE_DIR']
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks')
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')

INDIVS = glob_wildcards(os.path.join(DIR_FIRST_PASS, '{indiv}-clean.bam'))
INDIVS = INDIVS[0]
INDIVS.sort()
INDIVS = INDIVS[:10]

SUFFIX = 'combined_peaks_1_10.bam'

''' Snakemake rules '''
rule all:
    input:
        os.path.join(DIR_FIRST_PASS, SUFFIX)

rule merge_bam:
    input:
        expand(os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam'), indiv = INDIVS)
    output:
        os.path.join(DIR_FIRST_PASS, SUFFIX)
    shell:
        'samtools merge {output} {input}'

rule call_peaks:
    input:
        os.path.join(DIR_FIRST_PASS, SUFFIX)
    output:
        os.path.join(PEAK_DIR, 'combined_peaks.narrowPeak')
    params:
        outdir =  PEAK_DIR,
        prefix =  'combined_peaks',
    conda:
        "envs/env_py37.yml"
    shell:
        """
        {MACS2} callpeak -t {input} -f BAM -g hs --shift -75 --extsize 150 --outdir {params.outdir} -n {params.prefix} -q 0.05
        """


