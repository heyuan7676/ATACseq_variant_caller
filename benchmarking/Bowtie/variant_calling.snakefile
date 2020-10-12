import os

## TODO: interesting_region

DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')
GENOME_DICT = '.'.join(GENOME_STAR.split('.')[:-1]) + '.dict'

'''
Variant calling using GATK
'''

rule GATK_haplotypecaller:
    input:
        bam = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-RG-dedup-cleanH.bam'),
        bai = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-RG-dedup-cleanH.bam.bai'),
        genome_dict = GENOME_DICT,
        interesting_region = os.path.join(DIR_FIRST_PASS, '{indiv}' + '_peaks.narrowPeak.bed') 
    output:
        gvcf_gz = os.path.join(BOWTIE_DIR, '{indiv}-RG-dedup-hapcal.g.vcf.gz')
    params:
        tmp = os.path.join(BOWTIE_DIR, 'tmp')
    threads: THREADS
    shell:
        '{GATK} --java-options "-XX:ParallelGCThreads={threads}" HaplotypeCaller -R {GENOME_STAR} -L {input.interesting_region} -ip 1000 \
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
        {BCFTOOLS} convert --gvcf2vcf {input} -f {GENOME_STAR} > {output}
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



