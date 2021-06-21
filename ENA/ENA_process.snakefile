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
    conda:
        "envs/env_py37.yml"
    shell:
        """
        cd {params.outdir}
        kingfisher get -r {wildcards.indiv} -m ena-ascp ena-ftp
        gzip {params.reads1}
        gzip {params.reads2} 
        """


rule align_to_hg38:
    input:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + '_1.fastq.gz'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + '_2.fastq.gz'),
    output:
        bam = temp(os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam')),
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
BOWTIE_FA = '/work-zfs/abattle4/heyuan/database/bowtie2/grch38/GRCh38_noalt_as/GRCh38_noalt_as.fa'
'''Call peaks using MACS2'''
rule peak_calling_macs2:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam'),
    output:
        bed = os.path.join(PEAK_DIR_MACS2, '{indiv}' + '-clean.bed'),
        output = os.path.join(PEAK_DIR_MACS2, '{indiv}' + '_peaks.narrowPeak'),
        cram = os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.cram')
    params:
        outdir =  PEAK_DIR_MACS2,
        prefix =  '{indiv}',
    conda:
        "envs/env_py37.yml"
    shell:
        """
        {BEDTOOLS} bamtobed -i {input} > {output.bed}
        {MACS2} callpeak -t {output.bed} -f BED -g hs --nomodel --shift -100 --extsize 200 --outdir {params.outdir} -n {params.prefix} -q 0.05
        {CRAMTOOLS} cram --input-bam-file {input} --reference-fasta-file {BOWTIE_FA} --output-cram-file {output.cram}
        """


