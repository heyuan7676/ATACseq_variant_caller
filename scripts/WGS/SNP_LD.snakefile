import os
configfile: 'config.yaml'

'''COmpute R2 from VCF files for the 1000 Genome Project'''

BOWTIE_DIR = config['BOWTIE_DIR']
GT_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/onek_genome_data'
OneK_GENOME_DIR = config['OneK_GENOME_DIR']

CHROM = config['CHROM']
VCFTOOLS = config['VCFTOOLS']
PLINK = config['PLINK']

rule all:
    input:
        expand(os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2.ld'), chromosome = CHROM),
        expand(os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2.ld.variantPairs.txt'), chromosome = CHROM)


rule compute_r2:
    input:
        os.path.join(GT_dir, 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.recode.vcf')
    output:
        os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2.ld')
    params:
        binary = os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.plink'),
        r2 = os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2')
    shell:
        """
        {PLINK} --vcf {input} --recode --out {params.binary}
        {PLINK} --file {params.binary} --r2 inter-chr --ld-window-r2 0.2 --out {params.r2}
        rm {params.binary}*
        """


rule format:
    input:
        os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2.ld')
    output:
        os.path.join(GT_dir, 'plink', 'ALL.chr' + '{chromosome}' + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.r2.ld.variantPairs.txt')
    shell:
        """
        awk "{{print \$1"'"_"'"\$2, \$4"'"_"'"\$5, \$7}}" {input} > {output}
        """

