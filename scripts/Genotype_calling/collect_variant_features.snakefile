import os



'''
Process VCF file to obtain INFO
'''

rule format_info:
    input:
        os.path.join(VCF_DIR, '{indiv}.recode.vcf.gz')
    output:
        os.path.join(VCF_DIR, 'minDP{minDP}', "{indiv}.filtered.recode.INFO.txt")
    shell:
        """
        {BCFTOOLS} query -f '%CHROM\t%POS\t[%DP\t%GQ]\n' {input} > {output}
        """



rule obtain_snp_bed_file:
    input:
        os.path.join(VCF_DIR, "minDP{minDP}", '{indiv}' + ".filtered.minDP{minDP}.recode.dosage_genotype.bed")
    output:
        os.path.join(VCF_DIR, "minDP{minDP}", '{indiv}' + ".filtered.minDP{minDP}.recode.variants.bed")
    shell:
        """
        paste <(cut -d'_' -f1 {input} ) <(cut -d'_' -f2 {input} | awk "{{print \$1 = \$1-1}}")  <(cut -d'_' -f2 {input} | awk "{{print \$1}}") <(awk "{{print \$1}}" {input} ) | sed "1d" | sort -k1,1 -k2,2n > {output}
        """


rule closest_variants_to_peak:
    input:
        variants = os.path.join(VCF_DIR, "minDP{minDP}", '{indiv}' + ".filtered.minDP{minDP}.recode.variants.bed"),
        peaks = os.path.join(PEAK_DIR_MACS2, "{indiv}_peaks.narrowPeak")
    params:
        os.path.join(PEAK_DIR_MACS2, "{indiv}_peaks.narrowPeak_sorted")
    output:
        os.path.join(VCF_DIR, "minDP{minDP}", '{indiv}' + ".filtered.minDP{minDP}.recode.variants.toMACS2Peaks.txt")
    shell:
        """
        sort -k1,1 -k2,2n {input.peaks} > {params}
        {BEDTOOLS} closest -a {input.variants} -b {params} -d | awk "{{print \$4, \$8, \$15}}"> {output}
        rm {params}
        """
