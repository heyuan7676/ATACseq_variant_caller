import os


'''
Variant calling using GATK
'''

oneK_variants_locations='/work-zfs/abattle4/heyuan/Variant_calling/datasets/onek_genome_data/1k_genome.maf005.variants.locations.bed'
rule GATK_haplotypecaller:
    input:
        bam = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam'),
        bai = os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam.bai'),
        interesting_region = {oneK_variants_locations}
    output:
        gvcf_gz = temp(os.path.join(VCF_DIR, '{indiv}.gvcf.gz'))
    threads: THREADS
    shell:
        """
        {GATK} HaplotypeCaller -R {GENOME} -L {input.interesting_region} \
            -I {input.bam} -O {output.gvcf_gz} --tmp-dir {TMP_DIR} --native-pair-hmm-threads {threads} \
            -ERC GVCF
        """


'''
Prepare vcf files for 1). collecting variant called; 2). use as input for imputation
'''
rule filter_no_reads:
    input:
        os.path.join(VCF_DIR, '{indiv}.gvcf.gz')
    params:
        prefix = os.path.join(VCF_DIR, '{indiv}.save'),
        fn = os.path.join(VCF_DIR, '{indiv}.save.recode.vcf')
    output:
        vcf=os.path.join(VCF_DIR, '{indiv}.save.recode.vcf.gz'),
    shell:
        """
        {VCFTOOLS} --gzvcf {input} --min-meanDP 2 --recode --recode-INFO-all --out {params.prefix}
        bgzip {params.fn}
        """



rule gvcf2vcf_gatk:
    input:
        os.path.join(VCF_DIR, '{indiv}.gvcf.gz')
    output:
        gzfile = os.path.join(VCF_DIR, '{indiv}.forimputation.vcf.gz')
    params:
        vcffile = os.path.join(VCF_DIR, '{indiv}.forimputation.vcf'),
        gvcf_index_file = os.path.join(VCF_DIR, '{indiv}.gvcf.gz.tbi')
    shell:
        """
        {GATK} GenotypeGVCFs -R {GENOME} -V {input} -O {params.vcffile}
        bgzip {params.vcffile}
        rm -f {params.gvcf_index_file}
        """


'''
Filter variants
'''

rule filterVCF_minDP:
    input:
        os.path.join(VCF_DIR, '{indiv}.save.recode.vcf.gz')
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
        bgzip {params.vcffile}
        """



