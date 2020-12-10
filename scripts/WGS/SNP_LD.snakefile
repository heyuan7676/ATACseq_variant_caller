import os
configfile: 'config.yaml'

'''COmpute R2 from VCF files for the 1000 Genome Project'''

BOWTIE_DIR = config['BOWTIE_DIR']
GT_dir = os.path.join(BOWTIE_DIR, 'Genotype')
OneK_GENOME_DIR = config['OneK_GENOME_DIR']

CHROM = config['CHROM']
CHROM = ['22']
VCFTOOLS = config['VCFTOOLS']
PLINK = config['PLINK']

rule all:
    input:
        #expand(os.path.join(GT_dir, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.vcf'), chromosome = CHROM)
        expand(os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2'), chromosome = CHROM)


rule extract_maf_005_variants:
    input:
        os.path.join(OneK_GENOME_DIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz')
    output:
        os.path.join(GT_dir, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.vcf')
    params:
        prefix = os.path.join(GT_dir, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005')
    shell:
        """
        {VCFTOOLS} --gzvcf {input} --maf 0.05 --recode --recode-INFO-all --out {params}
        """



rule compute_r2:
    input:
        os.path.join(GT_dir, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.vcf')
    output:
        os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2')
    params:
        binary = os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.plink')
    shell:
        """
        {PLINK} --vcf {input} --recode --out {params}
        {PLINK} --file {params} --r2 inter-chr --ld-window-r2 0.2 --out {output}
        """
