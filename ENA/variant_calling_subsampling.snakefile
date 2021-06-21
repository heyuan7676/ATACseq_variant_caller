import os

RGPL = 'illumina'
RGPU = 'Unknown'

'''Map paired-ended fastq reads to genome using bowtie2'''

'''Subset '''
FQ_subsampled_DIR = os.path.join(FQ_DIR, 'subsampling')
rule subset_fastq:
    input:
        reads1 = os.path.join(FQ_DIR, '{indiv}' + '_1.fastq.gz'),
        reads2 = os.path.join(FQ_DIR, '{indiv}' + '_2.fastq.gz'),
    output:
        reads1 = temp(os.path.join(FQ_subsampled_DIR, '{indiv}' + '_1.fastq.subset.gz')),
        reads2 = temp(os.path.join(FQ_subsampled_DIR, '{indiv}' + '_2.fastq.subset.gz')),
    params:
        tool_dir = '/work-zfs/abattle4/heyuan/tools/seqtk',
        reads1 = os.path.join(FQ_subsampled_DIR, '{indiv}' + '_1.fastq.subset'),
        reads2 = os.path.join(FQ_subsampled_DIR, '{indiv}' + '_2.fastq.subset'),
    shell:
        """
        cd {params.tool_dir}
        ./seqtk sample -s100 {input.reads1} 5000000 > {params.reads1}
        ./seqtk sample -s100 {input.reads2} 5000000 > {params.reads2}
        gzip {params.reads1}
        gzip {params.reads2}
        """

suffix = 'subset'
rule subset_fastq_alignment:
    input:
        reads1 = os.path.join(FQ_subsampled_DIR, '{indiv}' + '_1.fastq.subset.gz'),
        reads2 = os.path.join(FQ_subsampled_DIR, '{indiv}' + '_2.fastq.subset.gz'),
    output:
        sam = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '.sam')),
        sam_sorted = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix+ '.sorted.sam')),
        sambai = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix+ '.sorted.sam.bai')),
        bam = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix+ '.bam')),
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
        os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '.bam')
    output:
        bam = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '_RG.bam')),
        bai = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '_RG.bam.bai'))
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
        os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '.bam')
    output:
        rg_bam = temp(os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + '-RG.bam')),
        bam = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '_nochr.bam'),
        bai = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '_nochr.bam.bai')
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



BOWTIE_FA = '/work-zfs/abattle4/heyuan/database/bowtie2/grch38/GRCh38_noalt_as/GRCh38_noalt_as.fa'
rule octopus_fast:
    input:
        bam = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '_RG.bam'),
        bai = os.path.join(BOWTIE_DIR, 'subsampling', '{indiv}' + suffix + '_RG.bam.bai'),
    output:
         os.path.join(VCF_DIR, '{indiv}.octopus.fast.vcf')
    shell:
        """
        ~/miniconda3/bin/octopus --reference {BOWTIE_FA} --reads {input.bam} -o {output} --fast --threads
        """
