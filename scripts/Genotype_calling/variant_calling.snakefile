import os

GENOME_DICT = '.'.join(GENOME_STAR.split('.')[:-1]) + '.dict'

'''
Variant calling using GATK
'''

oneK_variants_locations='/work-zfs/abattle4/heyuan/Variant_calling/datasets/onek_genome_data/1k_genome.maf005.variants.locations.bed'
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



'''
Filter variants
'''
rule gvcf2vcf:
    input:
        os.path.join(VCF_DIR, '{indiv}.gvcf.gz')
    params:
        os.path.join(VCF_DIR, '{indiv}.vcf')
    output:
        vcf=os.path.join(VCF_DIR, '{indiv}.vcf.gz'),
    shell:
        """
        {BCFTOOLS} convert --gvcf2vcf {input} -f {GENOME_STAR} > {params}
        bgzip {params}
        """


rule filterVCF_minDP:
    input:
        os.path.join(VCF_DIR, '{indiv}.vcf.gz')
    output:
        os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.filtered.minDP' + '{minDP}' + '.recode.vcf.gz')
    params:
        minimumdp = '{minDP}',
        prefix = os.path.join(VCF_DIR,'minDP{minDP}','{indiv}' + '.filtered.minDP' + '{minDP}'),
        vcffile = os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.filtered.minDP' + '{minDP}' + '.recode.vcf'),
        logfile = os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.filtered.minDP' + '{minDP}' + '.log')
    shell:
        """
        {VCFTOOLS} --gzvcf {input} --min-meanDP {params.minimumdp} --recode --recode-INFO-all --out {params.prefix}
        rm {params.logfile}
        bgzip {params.vcffile}
        """



'''
Process VCF file to obtain INFO
'''

rule format_info:
    input:
        os.path.join(VCF_DIR, '{indiv}.vcf.gz')
    params:
        os.path.join(VCF_DIR, '{indiv}.filtered.recode.INFO.vcf')
    output:
        os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.formatted.vcf")
    shell:
        """
        {BCFTOOLS} query -f '%CHROM\t%POS\t[%DP\t%PL\t%GQ]\n' {input} > {params}
        paste -d"_" <(awk "{{print \$1}}" {params}) <(awk "{{print \$2,\$4}}" {params}) | sed "1d" | sort -k1,1 -k2,2 > {output}
        rm {params}
        """



