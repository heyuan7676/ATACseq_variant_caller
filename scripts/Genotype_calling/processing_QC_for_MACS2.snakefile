import os

GENOME_DICT = '.'.join(GENOME.split('.')[:-1]) + '.dict'
RGPL = 'illumina'
RGPU = 'Unknown'

rule build_genome_dict:
    '''
    Build `.dict` file for reference genome
    '''
    input:
        GENOME
    output:
        faidx = GENOME + '.fai',
        dict = GENOME_DICT
    shell:
        """
        {SAMTOOLS} faidx {input} -o {output.faidx}
        {PICARD} CreateSequenceDictionary R={input} O={output.dict}
        """


rule build_vcf_index:
    input:
        VCFFN
    output:
        VCFFN + '.idx'
    conda:
        "../envs/env_py37.yml"
    threads: THREADS
    shell:
        '{GATK} IndexFeatureFile -I {input}'



'''Remove reads tagged by picard Duplicate'''
rule removedup:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam')
    output:
        bam = temp(os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup.bam')),
        metric = os.path.join(BOWTIE_DIR, 'dedup_metrics', '{indiv}' + '.dedup.metrics'),
        bam_format = temp(os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup_nochr.bam')),
        bam_clean = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam'),
        bai_clean = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam.bai')
    threads: THREADS
    shell:
        """
        {PICARD} MarkDuplicates INPUT= {input} OUTPUT= {output.bam} METRICS_FILE= {output.metric} TMP_DIR={TMP_DIR}
        {SAMTOOLS} view -H {output.bam} | sed  -e 's/chr//g'  | {SAMTOOLS} reheader - {output.bam} > {output.bam_format}
        {SAMTOOLS} view -b -F 0x400 {output.bam_format} > {output.bam_clean}
        {SAMTOOLS} index {output.bam_clean} 
        """

