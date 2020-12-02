import os

configfile: 'config.yaml'

BOWTIE_DIR = config['BOWTIE_DIR']
VCF_DIR = os.path.join(BOWTIE_DIR, 'VCF_files')
IMPUTATION_DIR = os.path.join(BOWTIE_DIR, 'Imputation')

INDIVS = glob_wildcards(os.path.join(VCF_DIR, '{indiv}.gvcf.gz'))
INDIVS = INDIVS[0]
INDIVS = INDIVS[:2]

CHROM = config['CHROM']

wildcard_constraints:
    chr="\d+"

'''Collect genotype data'''

rule all:
    input:
        expand(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snpids"), indiv = INDIVS),
        os.path.join(IMPUTATION_DIR,  "union-SNPs_imputed.bed"),
        expand(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.filtered.genotype.txt_matrix"), indiv = INDIVS),
        expand(os.path.join(IMPUTATION_DIR, "gt_by_sample_matrix_chr{chr}.txt"), chr = CHROM)


rule merge_chrs:
    input:
        files = expand(os.path.join(IMPUTATION_DIR, "{{indiv}}", "chr{chr}.imputed.GRCh38.genotype.txt"), chr = CHROM)
    output:
        temp = temp(os.path.join(IMPUTATION_DIR, "{indiv}", "{indiv}.imputed.GRCh38.genotype.txt_temp")),
        snps = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snps"),
        genotype = os.path.join(IMPUTATION_DIR, "{indiv}", "{indiv}.imputed.GRCh38.genotype.txt"),
    shell:
        """
        cat {input.files} | sed 's/chr//g' > {output.temp}
        paste -d"_" <(awk "{{print \$1}}" {output.temp}) <(awk "{{print \$2}}" {output.temp}) > {output.snps}
        paste {output.snps} {output.temp} | sort -k1,1 -k2,2n > {output.genotype}
        awk "{{print \$1, \$5}}" {output.genotype} > {output.snps}
        """



rule keep_biallelic_variants:
    input:
        snps = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snps"),
    output:
        snp_ids = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snpids")
    shell:
        """
        awk "{{print \$1}}" {input.snps} | sort | uniq -u > {output.snp_ids}
        """



rule genotype_biallelic_variants:
    input:
        header = os.path.join(IMPUTATION_DIR, "{indiv}", "chr22.imputed.GRCh37.genotype.txt"),
        genotype = os.path.join(IMPUTATION_DIR, "{indiv}", "{indiv}.imputed.GRCh38.genotype.txt"),
        snp_ids = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snpids")
    params:
        fn = os.path.join(IMPUTATION_DIR, "{indiv}", "{indiv}.imputed.GRCh38.biallelic.genotype.txt_temp")
    output:
        gt = os.path.join(IMPUTATION_DIR, "{indiv}", "{indiv}.imputed.GRCh38.biallelic.genotype.txt"),
        snps = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.biallelic.genotype.txt_snps")
    shell:
        """
        awk -F '\t' "NR==FNR {{id[\$1]; next}} \$1 in id" {input.snp_ids} {input.genotype} | sort -k1,1 >> {params.fn}
        awk "{{print \$1,\$5}}" {params.fn} | sed 's/ /	/g' > {output.snps}
        head -n1 {input.header} > {output.gt}
        awk "{{print \$2,\$3,\$4,\$5}}" {params.fn} | sed 's/ /	/g' >> {output.gt}
        rm {params.fn}
        """


rule union_set_SNPs:
    input:
        expand(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.genotype.txt_snpids"), indiv = INDIVS)
    output:
        fn1 = os.path.join(IMPUTATION_DIR,  "union-SNPs_imputed.bed"),
        fn2 = os.path.join(IMPUTATION_DIR, 'union-SNPs_imputed.info.bed')
    shell:
        """
        sort -mu {input} > {output.fn1}
        paste <(cat {output.fn1}) <(cat {output.fn1} | cut -d'_' -f 1) <(cat {output.fn1} | cut -d'_' -f 2) > {output.fn2}
        """


rule collect_genotype_union:
    input:
        gt = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.imputed.GRCh38.biallelic.genotype.txt_snps"),
        snps = os.path.join(IMPUTATION_DIR,  "union-SNPs_imputed.bed")
    output:
        gt = os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.filtered.genotype.txt_matrix")
    shell:
        """
        join -e 0 -o auto -a 1 -j 1 {input.snps} {input.gt} | awk "{{print \$2}}" > {output.gt}
        """


rule obtain_matrix:
    input:
        gt = expand(os.path.join(IMPUTATION_DIR, "{indiv}","{indiv}.filtered.genotype.txt_matrix"), indiv = INDIVS),
        snps = os.path.join(IMPUTATION_DIR, 'union-SNPs_imputed.info.bed')
    output:
        bychr = os.path.join(IMPUTATION_DIR, "gt_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.gt} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """

