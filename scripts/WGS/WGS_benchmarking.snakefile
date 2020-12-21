import os
configfile: 'config.yaml'

'''Extract variants with MAF > 0.05 among samples from the 74 cohort'''

CHROM = config['CHROM']
VCFTOOLS = config['VCFTOOLS']

OneK_GENOME_DIR = config['OneK_GENOME_DIR']
GT_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'


rule all:
    input:
        expand(os.path.join(GT_DIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.vcf'),chromosome = CHROM)


SAMPLE_fn = 'samples_GBR.txt'
rule extract_maf_005_variants:
    input:
        os.path.join(OneK_GENOME_DIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz')
    output:
        os.path.join(GT_DIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.vcf'),
    params:
        prefix = os.path.join(GT_DIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005')
    shell:
        """
        {VCFTOOLS} --gzvcf {input} --keep {SAMPLE_fn} --maf 0.05 --recode --recode-INFO-all --out {params}
        """

