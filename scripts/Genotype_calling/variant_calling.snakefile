import os

GENOME_DICT = '.'.join(GENOME_STAR.split('.')[:-1]) + '.dict'

'''
Variant calling using GATK
'''

oneK_variants_locations='/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype/maf005/1k_genome.variants.locations.bed'
rule GATK_haplotypecaller:
    input:
        bam = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam'),
        bai = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam.bai'),
        genome_dict = GENOME_DICT,
        interesting_region = {oneK_variants_locations}
    output:
        gvcf_gz = temp(os.path.join(VCF_DIR, '{indiv}.gvcf.gz'))
    threads: THREADS
    shell:
        '{GATK} HaplotypeCaller -R {GENOME_STAR} -L {input.interesting_region} \
            -I {input.bam} -O {output.gvcf_gz} --tmp-dir {TMP_DIR} --native-pair-hmm-threads {threads} \
            -ERC GVCF'



rule gvcf2vcf_gatk:
    input:
        os.path.join(VCF_DIR, '{indiv}.gvcf.gz')
    output:
        vcffile = os.path.join(VCF_DIR, '{indiv}.vcf'),
        gzfile = os.path.join(VCF_DIR, '{indiv}.vcf.gz')
    shell:
        """
        {GATK} GenotypeGVCFs -R {GENOME_STAR} -V {input} -O {output.vcffile}
        bgzip {output.vcffile}
        """


'''
Process VCF file
'''

rule format_info:
    input:
        os.path.join(VCF_DIR, '{indiv}.vcf.gz')
    output:
        info=os.path.join(VCF_DIR, '{indiv}.filtered.recode.INFO.vcf'),
        fn1=temp(os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_t")),
        fn2=os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.formatted.vcf")
    shell:
        """
        {BCFTOOLS} query -f '%CHROM\t%POS\t[%DP\t%PL\t%GQ]\n' {input} > {output.info}
        awk "{{if((\$3>1)) print \$0}}" {input} > {output.fn1}
        paste -d"_" <(awk "{{print \$1}}" {output.fn1}) <(awk "{{print \$2,\$4}}" {output.fn1}) | sed "1d" | sort -k1,1 -k2,2 > {output.fn2}
        """


'''
Filter variants
'''


rule filterVCF_minDP:
    input:
        os.path.join(VCF_DIR, '{indiv}.vcf.gz')
    output:
        temp(os.path.join(VCF_DIR, '{indiv}' + '.filtered.minDP' + '{minDP}' + '.recode.vcf'))
    params:
        minimumdp = '{minDP}',
        prefix = os.path.join(VCF_DIR,'{indiv}' + '.filtered.minDP' + '{minDP}')
    shell:
        """
        {VCFTOOLS} --gzvcf {input} --min-meanDP {params.minimumdp} --recode --recode-INFO-all --out {params.prefix}
        """


'''
Call genotypes
'''

rule call_genotype:
    input:
        os.path.join(VCF_DIR, '{indiv}.filtered.minDP' + '{minDP}' + '.recode.vcf')
    output:
        os.path.join(GENOTYPE_DIR, 'minDP{minDP}', '{indiv}.filtered.genotype.minDP' + '{minDP}' + '.txt')
    conda:
        "../envs/env_py37.yml"
    shell:
        'vcf-to-tab < {input} > {output}'



