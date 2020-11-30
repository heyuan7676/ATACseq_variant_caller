import os

configfile: 'config.yaml'

BOWTIE_DIR = config['BOWTIE_DIR']
VCF_DIR = os.path.join(BOWTIE_DIR, 'VCF_files')
IMPUTATION_DIR = os.path.join(BOWTIE_DIR, 'Imputation')

INDIVS = glob_wildcards(os.path.join(VCF_DIR, '{indiv}.gvcf.gz'))
INDIVS = INDIVS[0]

CHROM = config['CHROM']

'''Collect genotype data'''

rule all:
    input:
        #expand(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snpids"), indiv = INDIVS)
        os.path.join(IMPUTATION_DIR,  "union-SNPs_imputed_v2.bed"),
        #expand(os.path.join(IMPUTATION_DIR, "gt_by_sample_matrix_chr{chr}.txt"), chr = CHROM)


rule collect_gt_each:
    input:
        os.path.join(IMPUTATION_DIR, "{indiv}", "{indiv}.imputed.GRCh38.genotype.txt")
    output:
        snps = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snps"),
        snp_ids = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snpids")
    shell:
        """
        paste -d"_" <(awk "{{print \$1}}" {input}) <(awk "{{print \$2,\$4}}" {input}) | sed "1d" | sort -k1,1 > {output.snps}
        awk "{{print \$1}}" {output.snps} > {output.snp_ids}
        """


rule union_set_SNPs:
    input:
        expand(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snpids"), indiv = INDIVS)
    output:
        fn1 = temp(os.path.join(IMPUTATION_DIR,  "union-SNPs_imputed.bed_temp_v2")),
        fn2 = os.path.join(IMPUTATION_DIR,  "union-SNPs_imputed_v2.bed")
    shell:
        """
        sort -mu {input} > {output.fn1}
        paste <( cat {output.fn1}) <(cat {output.fn1} | cut -d"_" -f1) <(cat {output.fn1}| cut -d"_" -f2) > {output.fn2}
        """


rule collect_genotype_union:
    input:
        gt = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snps"),
        snps = os.path.join(IMPUTATION_DIR,  "union-SNPs_imputed.bed")
    output:
        gt = temp(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.filtered.genotype.txt_matrix"))
    shell:
        """
        join -e0 -a 1 -a 2 -j 1 {input.snps} -o auto {input.gt} | awk "{{print \$4}}"> {output.gt}
        """



rule obtain_matrix:
    input:
        gt = expand(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.filtered.genotype.txt_matrix"), indiv = INDIVS),
        snps = os.path.join(IMPUTATION_DIR, 'union-SNPs_imputed.bed')
    output:
        bychr = os.path.join(IMPUTATION_DIR, "gt_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.gt} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """


