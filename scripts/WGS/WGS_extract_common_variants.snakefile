import os
configfile: 'config.yaml'

'''Extract variants with population MAF > 0.05 from the 1000 Genome Project'''

CHROM = config['CHROM']
VCFTOOLS = config['VCFTOOLS']

OneK_GENOME_DIR = config['OneK_GENOME_DIR']
OUTDIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/onek_genome_data'


rule all:
    input:
        expand(os.path.join(OUTDIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.vcf'), chromosome = CHROM),
        expand(os.path.join(OUTDIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.variants.txt'), chromosome = CHROM)


rule extract_maf_005_variants:
    input:
        os.path.join(OneK_GENOME_DIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz')
    output:
        os.path.join(OUTDIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.vcf'),
    params:
        prefix = os.path.join(OUTDIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005')
    shell:
        """
        {VCFTOOLS} --gzvcf {input} --maf 0.05 --recode --recode-INFO-all --out {params}
        """


rule extract_variants_ids:
    input:
        os.path.join(OUTDIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.vcf'),
    output:
        os.path.join(OUTDIR, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.variants.txt'),
    shell:
        """
        echo "CHR_POS #CHROM POS" > {output}
        awk "{{print \$1"'"_"'"\$2, \$1, \$2}}" {input} | grep -v "#" >> {output}
        """


