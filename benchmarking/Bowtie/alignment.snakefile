import os

'''Map paired-ended fastq reads to genome using bowtie2'''

rule align_to_hg38:
    input:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + '_r1.fixed.fastq.gz'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + '_r2.fixed.fastq.gz'),
    output:
        sam = temp(os.path.join(BOWTIE_DIR, '{indiv}' + '.bowtie2.grch38.sam'))
    params:
        index = BOWTIE_GENOME_INDEX
    threads:
        THREADS
    conda:
        "envs/env_py37.yml"
    shell:
        'bowtie2 -t -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam} --threads {threads}'


rule sort_sam:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + '.bowtie2.grch38.sam')
    output:
        temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam'))
    shell:
        '{SAMTOOLS} sort {input} > {output}'


rule index_sam:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam')
    output:
        temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam.bai'))
    shell:
        '{SAMTOOLS} index {input}'


rule keep_autosomal_genome:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam'),
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam.bai')
    output:
        temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam'))
    threads:
        THREADS
    shell:
        '{SAMTOOLS} view -@ {threads} -b {input} chr{{1..22}} > {output}'


