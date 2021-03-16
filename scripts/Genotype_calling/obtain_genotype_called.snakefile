
'''Collect genotype data'''

suffix = 'recode'
rule extract_snpids_called:
    input:
        genotype = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "{indiv}.filtered.minDP{minDP}." + suffix + ".dosage_genotype.bed")
    output:
        snp_ids = temp(os.path.join(VCF_DIR, 'minDP' + "{minDP}", "{indiv}.filtered.minDP{minDP}." + suffix + ".dosage_genotype.snpids.bed"))
    shell:
        """
        awk "{{print \$1}}" {input.genotype} | sort | uniq -u > {output.snp_ids}
        """


rule union_set_SNPs_called:
    input:
        expand(os.path.join(VCF_DIR, 'minDP' + "{{minDP}}", "{indiv}.filtered.minDP{{minDP}}." + suffix + ".dosage_genotype.snpids.bed"), indiv = INDIVS)
    output:
        fn1 = os.path.join(VCF_DIR, 'minDP' + "{minDP}",  "union-SNPs_called." + suffix + ".bed"),
        fn2 = os.path.join(VCF_DIR, 'minDP' + "{minDP}", 'union-SNPs_called.' + suffix + '.info.bed')
    shell:
        """
        sort -mu {input} > {output.fn1}
        paste <(cat {output.fn1}) <(cat {output.fn1} | cut -d'_' -f 1) <(cat {output.fn1} | cut -d'_' -f 2) > {output.fn2}
        """


rule collect_genotype_union_called:
    input:
        gt = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "{indiv}.filtered.minDP{minDP}." + suffix + ".dosage_genotype.bed"),
        snps = os.path.join(VCF_DIR, 'minDP' + "{minDP}",  "union-SNPs_called." + suffix + ".bed"),
    output:
        dosage = temp(os.path.join(VCF_DIR, 'minDP' + "{minDP}", "{indiv}." + suffix + ".dosage_genotype.matrix.txt")),
    shell:
        """
        join -e -1 -o auto -a 1 -j 1 {input.snps} <(sort -k1,1 {input.gt}) | awk "{{print \$3}}" > {output.dosage}
        """



rule obtain_dosage_matrix_called:
    input:
        dosage = expand(os.path.join(VCF_DIR, 'minDP' + "{{minDP}}", "{indiv}." + suffix + ".dosage_genotype.matrix.txt"), indiv = INDIVS),
        snps = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "union-SNPs_called." + suffix + ".info.bed")
    output:
        bychr = os.path.join(VCF_DIR, 'minDP' + "{minDP}", "dosage_by_sample_matrix_chr{chr}." + suffix + ".txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.dosage} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """


