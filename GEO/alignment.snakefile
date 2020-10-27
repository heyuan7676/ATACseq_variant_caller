import os

'''Map paired-ended fastq reads to genome using bowtie2'''

rule align_to_hg38:
    input:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + '_1.fastq'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + '_2.fastq'),
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
        sam=os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam'),
        sambai=os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam.bai')
    output:
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam')
    threads:
        THREADS
    shell:
        '{SAMTOOLS} view -@ {threads} -b {input.sam} chr{{1..22}} > {output}'


