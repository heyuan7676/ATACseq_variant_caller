
BLACKLIST = os.path.join(PEAK_DIR, 'hg38.blacklist.bed')

'''Call peaks using MACS2'''
rule peak_calling:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam')
    output:
        bed = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bed'),
        output = os.path.join(PEAK_DIR_MACS2, '{indiv}' + '_peaks.narrowPeak')
    params:
        outdir =  PEAK_DIR_MACS2,
        prefix =  '{indiv}',
    conda:
        "envs/env_py37.yml"
    shell:
        """
        {BEDTOOLS} bamtobed -i {input} > {output.bed}
        {MACS2} callpeak -t {output.bed} -f BED -g hs --nomodel --shift -100 --extsize 200 --max-gap 100 --outdir {params.outdir} -n {params.prefix} -q 0.05
        """


'''Call peaks using Genrich'''

rule sortBam:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam')
    output:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.sorted.bam')
    shell:
        """
        {SAMTOOLS} sort -n {input} > {output}
        """



rule peak_calling:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.sorted.bam')
    output:
        Signal=os.path.join(PEAK_Genrich_DIR, '{indiv}' + '_peaks.narrowPeak'),
        bedgraph=os.path.join(PEAK_Genrich_DIR, '{indiv}' + '_peaks.bedGraph')
    shell:
        """
        {GENRICH} -t {input} -j -r -o {output.Signal} -f {output.bedgraph} -v -E {BLACKLIST} -q 0.05
        """


