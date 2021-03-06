
rule collect_gt_each:
    input:
        os.path.join(GENOTYPE_DIR, "minDP{minDP}", "{indiv}.filtered.genotype.minDP{minDP}.txt")
    output:
        snps = os.path.join(GENOTYPE_DIR, "minDP{minDP}","{indiv}.filtered.genotype.txt_snps"),
    shell:
        """
        paste -d"_" <(awk "{{print \$1}}" {input}) <(awk "{{print \$2,\$4}}" {input}) | sed "1d" | sort -k1,1 > {output.snps}
        """


rule union_set_SNPs:
    input:
        expand(os.path.join(GENOTYPE_DIR, "minDP{{minDP}}","{indiv}.filtered.genotype.txt_snps"), indiv = INDIVS)
    output:
        fn1 = temp(os.path.join(GENOTYPE_DIR,  "minDP{minDP}", "union-SNPs_minDP{minDP}.bed_snps")),
        fn2 = os.path.join(GENOTYPE_DIR,  "minDP{minDP}", "union-SNPs_minDP{minDP}.bed")
    shell:
        """
        rm -f {output.fn1}
        rm -f {output.fn2}
        awk "{{print \$1}}" {input} >> {output.fn1}
        sort -i {output.fn1} | uniq > {output.fn2}

        mv {output.fn2} {output.fn1}
        paste <( cat {output.fn1}) <(cat {output.fn1} | cut -d"_" -f1) <(cat {output.fn1}| cut -d"_" -f2) > {output.fn2}
        """


rule collect_genotype_union:
    input:
        gt = os.path.join(GENOTYPE_DIR, "minDP{minDP}","{indiv}.filtered.genotype.txt_snps"),
        snps = os.path.join(GENOTYPE_DIR,  "minDP{minDP}", 'union-SNPs_minDP{minDP}.bed')
    output:
        gt = temp(os.path.join(GENOTYPE_DIR, "minDP{minDP}","{indiv}.filtered.genotype.txt_matrix"))
    shell:
        """
        join -e0 -a 1 -a 2 -j 1 {input.snps} -o auto {input.gt} | awk "{{print \$4}}"> {output.gt}
        """

rule obtain_matrix:
    input:
        gt = expand(os.path.join(GENOTYPE_DIR, "minDP{{minDP}}","{indiv}.filtered.genotype.txt_matrix"), indiv = INDIVS),
        snps = os.path.join(GENOTYPE_DIR,  "minDP{minDP}", 'union-SNPs_minDP{minDP}.bed')
    output:
        bychr = os.path.join(GENOTYPE_DIR, "minDP{minDP}", "gt_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.gt} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """



'''
Use the most comprehensive SNP set to collect INFO - avoid too many files
'''

rule collect_INFO_union:
    input:
        info = os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.formatted.vcf"),
        snps = os.path.join(GENOTYPE_DIR,  "minDP2", 'union-SNPs_minDP2.bed')
    output:
        info = os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_matrix")
    shell:
        """
        join -e0 -a 1 -a 2 -j 1 {input.snps} -o auto {input.info} | awk "{{print \$4}}"> {output.info}
        """

rule obtain_INFO_matrix:
    input:
        info = expand(os.path.join(VCF_DIR, "{indiv}.filtered.recode.INFO.vcf_matrix"), indiv = INDIVS),
        snps = os.path.join(GENOTYPE_DIR,  "minDP2", 'union-SNPs_minDP2.bed')
    output:
        bychr = os.path.join(VCF_DIR, "gt_info_by_sample_matrix_chr{chr}_atac.txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.info} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}" | sed 's/	/ /g' >> {output.bychr}
        """




