import os

RGPL = 'illumina'
RGPU = 'Unknown'

'''Map paired-ended fastq reads to genome using bowtie2'''

rule align_to_hg38_subset:
    input:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + 'subset_1.fastq.gz'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + 'subset_2.fastq.gz'),
    output:
        sam = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '.sam')),
        sam_sorted = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '.sorted.sam')),
        sambai = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '.sorted.sam.bai')),
        bam = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '.bam'),
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
        {SAMTOOLS} index {output.bam}
        """


rule add_RG:
    input:
        os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '.bam')
    output:
        bam = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '_RG.bam'),
        bai = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '_RG.bam.bai')
    params:
        label = '{indiv}',
        grid = '{indiv}',
        rgsm = '{indiv}'
    shell:
        """
        {PICARD} AddOrReplaceReadGroups I={input}  O={output.bam}  RGID={params.grid} RGLB={params.label} RGPL={RGPL} RGSM={params.rgsm} RGPU={RGPU} TMP_DIR={TMP_DIR}
        {SAMTOOLS} index {output.bam}
        """


rule remove_chrprefix:
    input:
        os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '.bam')
    output:
        rg_bam = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + '-RG.bam')),
        bam = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '_nochr.bam'),
        bai = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '_nochr.bam.bai')
    params:
        label = '{indiv}',
        grid = '{indiv}',
        rgsm = '{indiv}'
    shell:
        """
        {PICARD} AddOrReplaceReadGroups I={input}  O={output.rg_bam}  RGID={params.grid} RGLB={params.label} RGPL={RGPL} RGSM={params.rgsm} RGPU={RGPU} TMP_DIR={TMP_DIR}
        {SAMTOOLS} view -H {output.rg_bam} | sed  -e 's/chr//g'  | {SAMTOOLS} reheader - {output.rg_bam} > {output.bam}
        {SAMTOOLS} index {output.bam}
        """


rule octopus_fast:
    input:
        bam = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '_RG.bam'),
        bai = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + SUFFIX + '_RG.bam.bai'),
    output:
         os.path.join(VCF_DIR, '{indiv}.octopus.fast.vcf')
    shell:
        """
        ~/miniconda3/bin/octopus --reference {BOWTIE_FA} --reads {input.bam} -o {output} --fast --threads
        """

