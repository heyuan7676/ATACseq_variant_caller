import os

'''Map paired-ended fastq reads to genome using bowtie2'''

rule download:
    output:
        reads1 = temp(os.path.join(FQ_DIR, '{indiv}' + '_1.fastq.gz')),
        reads2 = temp(os.path.join(FQ_DIR, '{indiv}' + '_2.fastq.gz')),
    params:
        outdir = FQ_DIR,
        reads1 = temp(os.path.join(FQ_DIR, '{indiv}' + '_1.fastq')),
        reads2 = temp(os.path.join(FQ_DIR, '{indiv}' + '_2.fastq')),
    shell:
        """
        cd {params.outdir}
        kingfisher get -r {indiv} -m ena-ascp aws-http prefetch
        gzip {params.reads1}
        gzip {params.reads2} 
        """


rule align_to_hg38:
    input:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + '_1.fastq.gz'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + '_2.fastq.gz'),
    output:
        bam = os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam'),
        sam = temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sam')),
        sam_sorted = temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sorted.sam')),
        sambai = temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.sorted.sam.bai'))
    params:
        index = BOWTIE_GENOME_INDEX,
    threads:
        THREADS
    conda:
        "envs/env_py37.yml"
    shell:
        """
        bowtie2 -t -x {params.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam} --threads {threads}
        {SAMTOOLS} sort {output.sam} > {output.sam_sorted}
        {SAMTOOLS} index {output.sam_sorted}
        {SAMTOOLS} view -@ {threads} -b {output.sam_sorted} chr{{1..22}} > {output.bam}
        """


BLACKLIST = os.path.join('/work-zfs/abattle4/heyuan/Variant_calling/datasets/hg38.blacklist.bed')

'''Call peaks using MACS2'''
rule peak_calling_macs2:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam'),
    output:
        bed = os.path.join(PEAK_DIR_MACS2, '{indiv}' + '-clean.bed'),
        output = os.path.join(PEAK_DIR_MACS2, '{indiv}' + '_peaks.narrowPeak')
    params:
        outdir =  PEAK_DIR_MACS2,
        prefix =  '{indiv}',
    conda:
        "envs/env_py37.yml"
    shell:
        """
        {BEDTOOLS} bamtobed -i {input} > {output.bed}
        {MACS2} callpeak -t {output.bed} -f BED -g hs --nomodel --shift -100 --extsize 200 --outdir {params.outdir} -n {params.prefix} -q 0.05
        """



'''Subset '''
rule subset_fastq:
    input:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + '_1.fastq.gz'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + '_2.fastq.gz'),
    output:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + 'subset_1.fastq.gz'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + 'subset_2.fastq.gz'),
    params:
        tool_dir = '/work-zfs/abattle4/heyuan/tools/seqtk',
        reads1 = os.path.join(FQ_DIR, '{indiv}' + 'subset_1.fastq'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + 'subset_2.fastq'),
    shell:
        """
        cd {params.tool_dir}
        ./seqtk sample -s100 {input.reads1} 5000000 > {params.reads1}
        ./seqtk sample -s100 {input.reads2} 5000000 > {params.reads2}
        gzip {params.reads1}
        gzip {params.reads2}
        """
