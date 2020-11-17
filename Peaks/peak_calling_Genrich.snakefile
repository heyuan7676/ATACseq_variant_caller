import os

configfile: 'config.yaml'

'''Load from config file'''

SAMTOOLS = config['SAMTOOLS']
GENRICH = config['GENRICH']

BOWTIE_DIR = config['BOWTIE_DIR']
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks')
PEAK_Genrich_DIR = os.path.join(BOWTIE_DIR, 'Peaks_Genrich')

INDIVS = glob_wildcards(os.path.join(DIR_FIRST_PASS, '{indiv}-clean.bam'))
INDIVS = INDIVS[0]

BLACKLIST = os.path.join(PEAK_DIR, 'hg38.blacklist.bed')

''' Snakemake rules '''
rule all:
    input:
        expand(os.path.join(PEAK_Genrich_DIR, '{indiv}' + '_peaks.Genrich.narrowPeak'), indiv = INDIVS)


rule sortBam:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam')
    output:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.sorted.bam')
    shell:
        """
        {SAMTOOLS} sort -n {input} > {output}
        """



'''Call peaks using Genrich'''
rule peak_calling:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.sorted.bam')
    output:
        Signal=os.path.join(PEAK_Genrich_DIR, '{indiv}' + '_peaks.Genrich.narrowPeak'),
        bedgraph=os.path.join(PEAK_Genrich_DIR, '{indiv}' + '_peaks.Genrich.bedGraph')
    params:
        prefix =  '{indiv}',
    shell:
        """
        {GENRICH} -t {input} -j -r -o {output.Signal} -f {output.bedgraph} -v -E {BLACKLIST} -q 0.05
        """


