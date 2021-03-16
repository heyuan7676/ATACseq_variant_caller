
'''Collect genotype data'''


rule union_set_SNPs:
    input:
        expand(os.path.join(VCF_DIR, 'minDP' + "{{minDP}}", "{indiv}", "{indiv}.filtered.minDP" + "{{minDP}}" + ".imputed.commonVariants.dosage_genotype.snpids.bed"), indiv = INDIVS)
    output:
        fn1 = os.path.join(VCF_DIR, 'minDP' + "{minDP}",  "union-SNPs_imputed.bed"),
        fn2 = os.path.join(VCF_DIR, 'minDP' + "{minDP}", 'union-SNPs_imputed.info.bed')
    shell:
        """
        sort -mu {input} > {output.fn1}
        paste <(cat {output.fn1}) <(cat {output.fn1} | cut -d'_' -f 1) <(cat {output.fn1} | cut -d'_' -f 2) > {output.fn2}
        """


rule collect_genotype_union:
    input:
        gt = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "{indiv}.filtered.minDP{minDP}.imputed.dosage_genotype.bed"),
        snps = os.path.join(VCF_DIR, 'minDP' + "{minDP}",  "union-SNPs_imputed.bed")
    output:
        dosage = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "{indiv}","{indiv}.filtered.genotype.dosage_matrix.txt"),
        #gt = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "{indiv}","{indiv}.filtered.genotype.gt_matrix.txt"),
    shell:
        """
        join -e 0 -o auto -a 1 -j 1 {input.snps} {input.gt} | awk "{{print \$3}}" > {output.dosage}
        """


rule obtain_gt_matrix:
    '''not used'''
    input:
        gt = expand(os.path.join(VCF_DIR, 'minDP' + "{{minDP}}", "{indiv}","{indiv}.filtered.genotype.gt_matrix.txt"), indiv = INDIVS),
        snps = os.path.join(VCF_DIR, 'minDP' + "{minDP}", 'union-SNPs_imputed.info.bed')
    output:
        bychr = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "gt_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.gt} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """


rule obtain_dosage_matrix:
    input:
        dosage = expand(os.path.join(VCF_DIR, 'minDP' + "{{minDP}}", "{indiv}","{indiv}.filtered.genotype.dosage_matrix.txt"), indiv = INDIVS),
        snps = os.path.join(VCF_DIR, 'minDP' + "{minDP}", 'union-SNPs_imputed.info.bed')
    output:
        bychr = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "dosage_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.dosage} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """


