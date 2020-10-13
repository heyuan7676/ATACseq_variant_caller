import os

GENOME_DICT = '.'.join(GENOME_STAR.split('.')[:-1]) + '.dict'
RGPL = 'illumina'
RGPU = 'Unknown'

'''Peak calling'''

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
        {MACS2} callpeak -t {input} -f BAM -g hs --shift -75 --extsize 150 --outdir {params.outdir} -n {params.prefix} -q 0.05
        """


