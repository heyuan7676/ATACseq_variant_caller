import os

SUFFIX = '.bowtie2.grch38.sortedByCoord.out'
RGPL = 'illumina'
RGPU = 'Unknown'

rule add_rg:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam')
    output:
        os.path.join(BOWTIE_DIR, '{indiv}' + '-RG.bam')
    params:
        tmp = os.path.join(BOWTIE_DIR, 'tmp/indiv'),
        label = '{indiv}',
        grid = '{indiv}',
        rgsm = '{indiv}'
    shell:
        '{PICARD} AddOrReplaceReadGroups I={input}  O={output}  RGID={params.grid} RGLB={params.label} RGPL={RGPL} RGSM={params.rgsm} RGPU={RGPU} TMP_DIR={params.tmp}'

rule mark_dup:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + '-RG.bam')
    output:
        bam = os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup.bam'),
        metric = os.path.join(BOWTIE_DIR, 'dedup_metrics', '{indiv}' + '.dedup.metrics')
    params:
        tmp = os.path.join(BOWTIE_DIR, 'tmp')
    threads: THREADS
    shell:
        '{PICARD} MarkDuplicates INPUT= {input} OUTPUT= {output.bam} METRICS_FILE= {output.metric} TMP_DIR={params.tmp}'


rule clean_header:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup.bam')
    output:
        os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup-cleanH.bam')
    shell:
        """
        {SAMTOOLS} view -H {input} | grep -v chr[a-Z] | grep -v random > {output}.temp
        {SAMTOOLS} view {input} >> {output}.temp
        {SAMTOOLS} view -h -b -S {output}.temp > {output}
        rm {output}.temp
        """

rule index_post_removeDup:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup-cleanH.bam'),
    output:
        os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup-cleanH.bam.bai')
    shell:
        '{SAMTOOLS} index {input}'



'''
Variant calling using GATK
'''

GENOME_DICT = '.'.join(GENOME.split('.')[:-1]) + '.dict'
rule build_genome_dict:
    '''
    Build `.dict` file for reference genome
    '''
    input:
        GENOME
    output:
        GENOME_DICT
    shell:
        '{PICARD} CreateSequenceDictionary R={input} O={output}'

rule GATK_haplotypecaller:
    input:
        bam = os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup-cleanH.bam'),
        bai = os.path.join(BOWTIE_DIR, '{indiv}' + '-RG-dedup-cleanH.bam.bai'),
        genome_dict = GENOME_DICT
    output:
        gvcf_gz = temp(os.path.join(BOWTIE_DIR, '{indiv}-RG-dedup-hapcal.g.vcf.gz'))
    params:
        tmp = os.path.join(BOWTIE_DIR, 'tmp')
    threads: THREADS
    shell:
        '{GATK} --java-options "-XX:ParallelGCThreads={threads}" HaplotypeCaller -R {GENOME} \
            -I {input.bam} -O {output.gvcf_gz} --tmp-dir {params.tmp} --native-pair-hmm-threads {threads} \
            -ERC GVCF'

rule filtergVCF_DP2:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}-RG-dedup-hapcal.g.vcf.gz')
    output:
        temp(os.path.join(BOWTIE_DIR, '{indiv}-RG-dedup-hapcal.filtered.gvcf.recode.vcf'))
    params:
        prefix = os.path.join(BOWTIE_DIR, '{indiv}-RG-dedup-hapcal.filtered.gvcf')
    shell:
        '{VCFTOOLS} --gzvcf {input} --min-meanDP 2 --recode --recode-INFO-all --out {params.prefix}'



'''
Process VCF file
'''

rule gvcf2vcf:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}-RG-dedup-hapcal.filtered.gvcf.recode.vcf')
    output:
        temp(os.path.join(BOWTIE_DIR, '{indiv}-RG-dedup-hapcal.filtered.recode.vcf'))
    shell:
        """
        {BCFTOOLS} convert --gvcf2vcf {input} -f {GENOME} > {output}
        """


oneK_variants_locations='/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/1k_genome.variants.locations.txt'
rule filterVCF_1k_variants:
    input:
        os.path.join(BOWTIE_DIR, '{indiv}-RG-dedup-hapcal.filtered.recode.vcf')
    output:
        os.path.join(BOWTIE_DIR, 'VCF_files', '{indiv}-RG-dedup-hapcal.filtered.variants.recode.vcf'),
        os.path.join(BOWTIE_DIR, 'VCF_files', '{indiv}-RG-dedup-hapcal.filtered.variants.recode.INFO.txt')
    params:
        prefix = os.path.join(BOWTIE_DIR, 'VCF_files', '{indiv}-RG-dedup-hapcal.filtered.variants')
    shell:
        """
        {VCFTOOLS} --vcf {input} --out {params.prefix} --positions {oneK_variants_locations} --recode
        {BCFTOOLS} query -f '%CHROM\t%POS\t[%DP\t%PL\t%GQ]\n' {params.prefix}.recode.vcf > {params.prefix}.recode.INFO.txt
        """


rule filterVCF_minDP:
    input:
        os.path.join(BOWTIE_DIR, 'VCF_files', '{indiv}' + '-RG-dedup-hapcal.filtered.variants.recode.vcf')
    output:
        temp(os.path.join(BOWTIE_DIR, 'VCF_files', '{indiv}' + '-RG-dedup-hapcal.filtered.variants.minDP' + '{minDP}' + '.recode.vcf'))
    params:
        minimumdp = '{minDP}',
        prefix = os.path.join(BOWTIE_DIR, 'VCF_files','{indiv}' + '-RG-dedup-hapcal.filtered.variants.minDP' + '{minDP}')
    shell:
        """
        {VCFTOOLS} --vcf {input} --min-meanDP {params.minimumdp} --recode --recode-INFO-all --out {params.prefix}
        """


'''
Call genotypes
'''

rule call_genotype:
    input:
        os.path.join(BOWTIE_DIR, 'VCF_files', '{indiv}-RG-dedup-hapcal.filtered.variants.minDP' + '{minDP}' + '.recode.vcf')
    output:
        os.path.join(BOWTIE_DIR, 'Called_GT', '{indiv}-RG-dedup-hapcal.filtered.genotype.minDP' + '{minDP}' + '.txt')
    conda:
        "envs/env_py37.yml"
    shell:
        'vcf-to-tab < {input} > {output}'


