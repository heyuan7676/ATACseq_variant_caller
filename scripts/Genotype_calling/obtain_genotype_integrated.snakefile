
'''Collect genotype data'''

rule extract_snpids_intg:
    input:
        genotype = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", "{indiv}_minDP{minDP}_variants_" + SUFFIX + ".txt")
    params:
        os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", "{indiv}_minDP{minDP}_variants_" + SUFFIX + ".txt_sorted")
    output:
        snp_ids = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", "{indiv}_minDP{minDP}_variants_" + SUFFIX + ".snpids.txt")
    shell:
        """
        sort -k1,1 {input.genotype} > {params}
        mv {params} {input.genotype}
        awk "{{print \$1}}" {input.genotype} | sort | uniq -u > {output.snp_ids}
        """


rule union_set_SNPs_intg:
    input:
        expand(os.path.join(INTERGRATED_DIR, 'minDP' + "{{minDP}}", "{indiv}_minDP{{minDP}}_variants_" + SUFFIX + ".snpids.txt"), indiv = INDIVS)
    output:
        fn1 = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}",  "union-SNPs_integreated." + SUFFIX + ".bed"),
        fn2 = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", 'union-SNPs_integreated.' + SUFFIX + '.info.bed')
    shell:
        """
        sort -mu {input} > {output.fn1}
        paste <(cat {output.fn1}) <(cat {output.fn1} | cut -d'_' -f 1) <(cat {output.fn1} | cut -d'_' -f 2) > {output.fn2}
        """


rule collect_genotype_union_intg:
    input:
        gt = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", "{indiv}_minDP{minDP}_variants_" + SUFFIX + ".txt"),
        snps = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}",  "union-SNPs_integreated." + SUFFIX + ".bed"),
    output:
        dosage = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", "{indiv}." + SUFFIX + ".dosage_genotype.matrix.txt"),
    shell:
        """
        join -e -1 -o auto -a 1 -j 1 {input.snps}  {input.gt} | awk "{{print \$2}}" > {output.dosage}
        """



rule obtain_dosage_matrix_intg:
    input:
        dosage = expand(os.path.join(INTERGRATED_DIR, 'minDP' + "{{minDP}}", "{indiv}." + SUFFIX + ".dosage_genotype.matrix.txt"), indiv = INDIVS),
        snps = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", 'union-SNPs_integreated.' + SUFFIX + '.info.bed')
    output:
        bychr = os.path.join(INTERGRATED_DIR, 'minDP' + "{minDP}", "dosage_by_sample_matrix_chr{chr}." + SUFFIX + ".txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.dosage} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """


