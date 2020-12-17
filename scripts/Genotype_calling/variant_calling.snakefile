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
        gzfile = os.path.join(VCF_DIR, '{indiv}.vcf.gz')
    params:
        vcffile = os.path.join(VCF_DIR, '{indiv}.vcf'),
        gvcf_index_file = os.path.join(VCF_DIR, '{indiv}.gvcf.gz.tbi')
    shell:
        """
        {GATK} GenotypeGVCFs -R {GENOME_STAR} -V {input} -O {params.vcffile}
        bgzip {params.vcffile}
        rm {params.gvcf_index_file}
        """


'''
Process VCF file
'''

rule format_info:
    input:
        os.path.join(VCF_DIR, '{indiv}.vcf.gz')
    output:
        info=os.path.join(VCF_DIR, '{indiv}.filtered.recode.INFO.vcf'),
        filename=os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.formatted.vcf")
    shell:
        """
        {BCFTOOLS} query -f '%CHROM\t%POS\t[%DP\t%PL\t%GQ]\n' {input} > {output.info}
        paste -d"_" <(awk "{{print \$1}}" {output.info}) <(awk "{{print \$2,\$4}}" {output.info}) | sed "1d" | sort -k1,1 -k2,2 > {output.filename}
        """


'''
Filter variants
'''


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
Call genotypes
'''

rule call_genotype:
    input:
        os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.filtered.minDP' + '{minDP}' + '.recode.vcf.gz')
    params:
        os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.filtered.minDP' + '{minDP}' + '.recode.vcf')
    output:
        os.path.join(GENOTYPE_DIR, 'minDP{minDP}', '{indiv}.filtered.genotype.minDP' + '{minDP}' + '.txt')
    conda:
        "../envs/env_py37.yml"
    shell:
        """
        gunzip {input}
        vcf-to-tab < {params} > {output}
        bgzip {params}
        """




