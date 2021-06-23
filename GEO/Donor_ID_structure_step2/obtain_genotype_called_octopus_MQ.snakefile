
'''Collect genotype data'''

suffix = 'recode'
rule MQ_octopus_extract_snpids_called:
    input:
        genotype = os.path.join(VCF_DIR, 'minMQ' + "{minMQ}", "{indiv}.filtered.minMQ{minMQ}." + suffix + ".dosage_genotype.octopus.bed")
    output:
        snp_ids = temp(os.path.join(VCF_DIR, 'minMQ' + "{minMQ}", "{indiv}.filtered.minMQ{minMQ}." + suffix + ".dosage_genotype.snpids.octopus.bed"))
    shell:
        """
        awk "{{print \$1}}" {input.genotype} | sort | uniq -u > {output.snp_ids}
        """


rule MQ_octopus_union_set_SNPs_called:
    input:
        expand(os.path.join(VCF_DIR, 'minMQ' + "{{minMQ}}", "{indiv}.filtered.minMQ{{minMQ}}." + suffix + ".dosage_genotype.snpids.octopus.bed"), indiv = INDIVS)
    output:
        fn1 = os.path.join(VCF_DIR, 'minMQ' + "{minMQ}",  "union-SNPs_called." + suffix + ".octopus.bed"),
        fn2 = os.path.join(VCF_DIR, 'minMQ' + "{minMQ}", 'union-SNPs_called.' + suffix + '.info.octopus.bed')
    shell:
        """
        sort -mu {input} > {output.fn1}
        paste <(cat {output.fn1}) <(cat {output.fn1} | cut -d'_' -f 1) <(cat {output.fn1} | cut -d'_' -f 2) > {output.fn2}
        """


rule MQ_octopus_collect_genotype_union_called:
    input:
        gt = os.path.join(VCF_DIR, 'minMQ' + "{minMQ}", "{indiv}.filtered.minMQ{minMQ}." + suffix + ".dosage_genotype.octopus.bed"),
        snps = os.path.join(VCF_DIR, 'minMQ' + "{minMQ}",  "union-SNPs_called." + suffix + ".octopus.bed"),
    output:
        dosage = os.path.join(VCF_DIR, 'minMQ' + "{minMQ}", "{indiv}." + suffix + ".dosage_genotype.matrix.txt"),
    shell:
        """
        join -e -1 -o auto -a 1 -j 1 {input.snps} <(sort -k1,1 {input.gt}) | awk "{{print \$3}}" > {output.dosage}
        """



rule MQ_octopus_obtain_dosage_matrix_called:
    input:
        dosage = expand(os.path.join(VCF_DIR, 'minMQ' + "{{minMQ}}", "{indiv}." + suffix + ".dosage_genotype.matrix.txt"), indiv = INDIVS),
        snps = os.path.join(VCF_DIR, 'minMQ' + "{minMQ}", "union-SNPs_called." + suffix + ".info.octopus.bed")
    output:
        bychr = os.path.join(VCF_DIR, 'minMQ' + "{minMQ}", "dosage_by_sample_matrix_chr{chr}." + suffix + ".txt")
    shell:
        """
        echo "CHR_POS CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.snps} {input.dosage} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}"  | sed 's/	/ /g' >> {output.bychr}
        """


