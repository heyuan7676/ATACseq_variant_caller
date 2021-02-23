import os

configfile: 'config.yaml'

'''Load from config file'''

BCFTOOLS = config['BCFTOOLS']
VCFTOOLS = config['VCFTOOLS']
BEDTOOLS = config['BEDTOOLS']

INDIVS = []
fn = open('Gencove_samples.txt', 'r')
for line in fn.readlines():
    INDIVS.append(line.rstrip())


GENCOVE_DIR = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Gencove'
oneK_variants = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/onek_genome_data/1k_genome.maf005.variants.locations.bed'

rule all:
    input:
        expand(os.path.join(GENCOVE_DIR, '{indiv}.gencove.bed'), indiv = INDIVS)


rule obtain_genotype_dosage:
    input:
        os.path.join(GENCOVE_DIR, '{indiv}.vcf.gz')
    output:
        mid = os.path.join(GENCOVE_DIR, '{indiv}.bed.mid'),
        dosage = os.path.join(GENCOVE_DIR, '{indiv}.gencove.bed')
    conda:
        "../scripts/envs/env_py37.yml"
    shell:
        """
        {BCFTOOLS} query -f '%CHROM\t%POS\t[%GT\t%DS]\n' {input} | awk "{{print \$1,\$2,\$4}}" | awk "{{print \$1"'"_"'"\$2, \$3}}" | sed 's/chr//g' > {output.mid}
        awk 'NR==FNR{{inFileA[$4]; next}} ($1 in inFileA)' {oneK_variants} {output.mid} > {output.dosage}
        """



