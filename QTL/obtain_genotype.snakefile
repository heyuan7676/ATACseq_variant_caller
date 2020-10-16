import os

configfile: 'config.yaml'

BOWTIE_DIR = config['BOWTIE_DIR']
GENOTYPE_DIR = os.path.join(BOWTIE_DIR, 'Called_GT')
VCF_DIR = os.path.join(BOWTIE_DIR, 'VCF_files')

INDIVS = glob_wildcards(os.path.join(GENOTYPE_DIR, '{indiv}.filtered.genotype.minDP3.txt'))
INDIVS = INDIVS[0]
CHROM = config['CHROM']

'''Collect genotype data'''

UNION_fn = 'union-SNPs.bed'
GT_fn = 'gt_by_sample_matrix.txt'

rule all:
    input:
        expand(os.path.join(GENOTYPE_DIR, "{indiv}.filtered.genotype.minDP3.txt_temp"), indiv = INDIVS),
        os.path.join(GENOTYPE_DIR, UNION_fn),
        expand(os.path.join(GENOTYPE_DIR, "{indiv}.filtered.genotype.minDP3.txt_matrix"), indiv = INDIVS),
        expand(os.path.join(GENOTYPE_DIR, "gt_by_sample_matrix_chr{chr}.txt"), chr = CHROM),
        expand(os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_temp"), indiv = INDIVS),
        expand(os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_matrix"), indiv = INDIVS),
        expand(os.path.join(VCF_DIR, "gt_info_by_sample_matrix_chr{chr}.txt"), chr = CHROM)


rule collect_gt_each:
    input:
        os.path.join(GENOTYPE_DIR, "{indiv}.filtered.genotype.minDP3.txt")
    output:
        snps = os.path.join(GENOTYPE_DIR, "{indiv}.filtered.genotype.minDP3.txt_temp"),
    shell:
        """
        paste -d"_" <(awk "{{print \$1}}" {input}) <(awk "{{print \$2,\$4}}" {input}) | sed "1d" | sort -k1,1 > {output.snps}
        """


rule union_set_SNPs:
    input:
        expand(os.path.join(GENOTYPE_DIR, "{indiv}.filtered.genotype.minDP3.txt_temp"), indiv = INDIVS)
    output:
        os.path.join(GENOTYPE_DIR, UNION_fn)
    params:
        os.path.join(GENOTYPE_DIR, UNION_fn + '_temp')
    shell:
        """
        rm -f {output}
        rm -f {params}
        awk "{{print \$1}}" {input} >> {params}
        sort -i {params} | uniq > {output}

        mv {output} {params}
        paste <( cat {params}) <(cat {params} | cut -d"_" -f1) <(cat {params}| cut -d"_" -f2) > {output}
        """


rule collect_genotype_union:
    input:
        gt = os.path.join(GENOTYPE_DIR, "{indiv}.filtered.genotype.minDP3.txt_temp"),
        snps = os.path.join(GENOTYPE_DIR, UNION_fn)
    output:
        gt = os.path.join(GENOTYPE_DIR, "{indiv}.filtered.genotype.minDP3.txt_matrix")
    shell:
        """
        join -e0 -a 1 -a 2 -j 1 {input.snps} -o auto {input.gt} | awk "{{print \$4}}"> {output.gt}
        """

rule obtain_matrix:
    input:
        gt = expand(os.path.join(GENOTYPE_DIR, "{indiv}.filtered.genotype.minDP3.txt_matrix"),indiv = INDIVS),
        snps = os.path.join(GENOTYPE_DIR, UNION_fn)
    output:
        bychr = os.path.join(GENOTYPE_DIR, "gt_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        paste {input.snps} {input.gt} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """


rule obtain_info:
    input:
        os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf")
    params:
        os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_t")
    output:
        os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_temp")
    shell:
        """
        awk "{{if((\$3>2)) print \$0}}" {input} > {params}
        paste -d"_" <(awk "{{print \$1}}" {params}) <(awk "{{print \$2,\$4}}" {params}) | sed "1d" | sort -k1,1 -k2,2 > {output}
        """


rule collect_INFO_union:
    input:
        info = os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_temp"),
        snps = os.path.join(GENOTYPE_DIR, UNION_fn)
    output:
        info = os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_matrix")
    shell:
        """
        join -e0 -a 1 -a 2 -j 1 {input.snps} -o auto {input.info} | awk "{{print \$4}}"> {output.info}
        """

rule obtain_INFO_matrix:
    input:
        info = expand(os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_matrix"), indiv = INDIVS),
        snps = os.path.join(GENOTYPE_DIR, UNION_fn)
    output:
        bychr = os.path.join(VCF_DIR, "gt_info_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        paste {input.snps} {input.info} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}" | sed 's/	/ /g' >> {output.bychr}
        """






