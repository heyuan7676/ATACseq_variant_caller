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


rule add_rg:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam')
    output:
        temp(os.path.join(BOWTIE_DIR, '{indiv}' + '-RG.bam'))
    params:
        label = '{indiv}',
        grid = '{indiv}',
        rgsm = '{indiv}'
    shell:
        '{PICARD} AddOrReplaceReadGroups I={input}  O={output}  RGID={params.grid} RGLB={params.label} RGPL={RGPL} RGSM={params.rgsm} RGPU={RGPU} TMP_DIR={TMP_DIR}'

rule mark_dup:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + '-RG.bam')
    output:
        bam = temp(os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup.bam')),
        metric = os.path.join(BOWTIE_DIR, 'dedup_metrics', '{indiv}' + '.dedup.metrics')
    threads: THREADS
    shell:
        '{PICARD} MarkDuplicates INPUT= {input} OUTPUT= {output.bam} METRICS_FILE= {output.metric} TMP_DIR={TMP_DIR}'


rule remove_chrprefix:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup.bam')
    output:
        temp(os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup_nochr.bam'))
    shell:
        """
        {SAMTOOLS} view -H {input} | sed  -e 's/chr//g'  | {SAMTOOLS} reheader - {input} > {output}
        """


rule build_bqsr_table:
    input:
        bam = os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup_nochr.bam'),
        genome = GENOME,
        known_vcf = VCFFN,
        genome_dict = GENOME_DICT,
        vcf_index = VCFFN + '.idx'
    output:
        table = temp(os.path.join(DIR_FIRST_PASS, '{indiv}' + '-RG-dedup.bqsr.table'))
    conda:
        "../envs/env_py37.yml"
    threads: THREADS
    shell:
        '{GATK} BaseRecalibrator -R {input.genome} -I {input.bam} --known-sites {input.known_vcf} -O {output.table}'


rule apply_bqsr:
    input:
        bam = os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup_nochr.bam'),
        genome = GENOME,
        table = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-RG-dedup.bqsr.table'),
    output:
        bam = temp(os.path.join(DIR_FIRST_PASS, '{indiv}' + '-RG-dedup-bqsr.bam')),
    threads: THREADS
    shell:
        '{GATK} ApplyBQSR -R {input.genome} -I {input.bam} --bqsr-recal-file {input.table} -O {output.bam}'



'''Remove reads tagged by picard Duplicate'''
rule removedup:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-RG-dedup-bqsr.bam')
    output:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam')
    shell:
        '{SAMTOOLS} view -b -F 0x400 {input} > {output}'


rule index_post_removeDup:
    input:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam'),
    output:
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam.bai')
    shell:
        '{SAMTOOLS} index {input}'

