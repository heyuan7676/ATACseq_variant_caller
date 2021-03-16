import os

'''Map paired-ended fastq reads to genome using bowtie2'''

rule align_to_hg38:
    input:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + '_1.fastq.gz'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + '_2.fastq.gz'),
    output:
        temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam'))
    params:
        index = BOWTIE_GENOME_INDEX,
        sam = temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam')),
        sam_sorted = temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sorted.sam')),
        sambai = temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sorted.sam.bai'))
    threads:
        THREADS
    conda:
        "envs/env_py37.yml"
    shell:
        """
        bowtie2 -t -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {params.sam} --threads {threads}
        {SAMTOOLS} sort {params.sam} > {params.sam_sorted}
        {SAMTOOLS} index {params.sam_sorted}
        {SAMTOOLS} view -@ {threads} -b {params.sam_sorted} chr{{1..22}} > {output}
        """
