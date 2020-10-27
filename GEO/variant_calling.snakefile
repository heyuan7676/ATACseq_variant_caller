import os

## TODO: interesting_region

GENOME_DICT = '.'.join(GENOME_STAR.split('.')[:-1]) + '.dict'

'''
Variant calling using GATK
'''

oneK_variants_locations='/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/1k_genome.variants.locations.bed'
#oneK_variants_locations='/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/temp.bed'
rule GATK_haplotypecaller:
    input:
        bam = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam'),
        bai = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam.bai'),
        genome_dict = GENOME_DICT,
        interesting_region = {oneK_variants_locations}
    output:
        gvcf_gz = os.path.join(VCF_DIR, '{indiv}.gvcf.gz')
    threads: THREADS
    shell:
        '{GATK} HaplotypeCaller -R {GENOME_STAR} -L {input.interesting_region} \
            -I {input.bam} -O {output.gvcf_gz} --tmp-dir {TMP_DIR} --native-pair-hmm-threads {threads} \
            -ERC GVCF'



rule filtergVCF_DP2:
    input:
        os.path.join(VCF_DIR, '{indiv}.gvcf.gz')
    output:
        temp(os.path.join(VCF_DIR, '{indiv}.filtered.gvcf.recode.vcf'))
    params:
        prefix = os.path.join(VCF_DIR, '{indiv}.filtered.gvcf')
    shell:
        '{VCFTOOLS} --gzvcf {input} --min-meanDP 2 --recode --recode-INFO-all --out {params.prefix}'



'''
Process VCF file
'''

rule gvcf2vcf:
    input:
        os.path.join(VCF_DIR, '{indiv}.filtered.gvcf.recode.vcf')
    output:
        vcf=os.path.join(VCF_DIR, '{indiv}.filtered.recode.vcf'),
        info=os.path.join(VCF_DIR, '{indiv}.filtered.recode.INFO.vcf')
    shell:
        """
        {BCFTOOLS} convert --gvcf2vcf {input} -f {GENOME_STAR} > {output.vcf}
        {BCFTOOLS} query -f '%CHROM\t%POS\t[%DP\t%PL\t%GQ]\n' {output.vcf} > {output.info}
        """


rule filterVCF_minDP:
    input:
        os.path.join(VCF_DIR ,'{indiv}' + '.filtered.recode.vcf')
    output:
        temp(os.path.join(VCF_DIR, '{indiv}' + '.filtered.minDP' + '{minDP}' + '.recode.vcf'))
    params:
        minimumdp = '{minDP}',
        prefix = os.path.join(VCF_DIR,'{indiv}' + '.filtered.minDP' + '{minDP}')
    shell:
        """
        {VCFTOOLS} --vcf {input} --min-meanDP {params.minimumdp} --recode --recode-INFO-all --out {params.prefix}
        """

rule filterVCF_GQ:
    input:
        os.path.join(VCF_DIR ,'{indiv}' + '.filtered.recode.vcf')
    output:
        temp(os.path.join(VCF_DIR, '{indiv}' + '.filtered.GQ' + '{GQ}' + '.recode.vcf'))
    params:
        minimumgq = '{GQ}',
        prefix = os.path.join(VCF_DIR,'{indiv}' + '.filtered.GQ' + '{GQ}')
    shell:
        """
        {VCFTOOLS} --vcf {input} --minGQ {params.minimumgq} --recode --recode-INFO-all --out {params.prefix}
        """



'''
Call genotypes
'''

rule call_genotype:
    input:
        os.path.join(VCF_DIR, '{indiv}.filtered.minDP' + '{minDP}' + '.recode.vcf')
    output:
        os.path.join(GENOTYPE_DIR, '{indiv}.filtered.genotype.minDP' + '{minDP}' + '.txt')
    shell:
        'vcf-to-tab < {input} > {output}'


rule call_genotype_GQ:
    input:
        os.path.join(VCF_DIR, '{indiv}.filtered.GQ' + '{GQ}' + '.recode.vcf')
    output:
        os.path.join(GENOTYPE_DIR, '{indiv}.filtered.genotype.GQ' + '{GQ}' + '.txt')
    shell:
        'vcf-to-tab < {input} > {output}'

