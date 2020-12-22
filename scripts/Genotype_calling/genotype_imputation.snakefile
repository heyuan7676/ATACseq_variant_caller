CHAIN_FN = "/work-zfs/abattle4/heyuan/tools/liftOver/hg38ToHg19.over.chain.gz"
REF_FN = "/work-zfs/abattle4/heyuan/database/GRCh37_reference/hg19.fa"


'''
Filter variants and prepare data for imputation
'''

rule gvcf2vcf_gatk:
    input:
        os.path.join(VCF_DIR, '{indiv}.gvcf.gz')
    output:
        gzfile = temp(os.path.join(VCF_DIR, '{indiv}.forimputation.vcf.gz'))
    params:
        vcffile = os.path.join(VCF_DIR, '{indiv}.forimputation.vcf'),
        gvcf_index_file = os.path.join(VCF_DIR, '{indiv}.gvcf.gz.tbi')
    shell:
        """
        {GATK} GenotypeGVCFs -R {GENOME_STAR} -V {input} -O {params.vcffile}
        bgzip {params.vcffile}
        rm {params.gvcf_index_file}
        """


rule filterVCF_minDP_for_imputation:
    input:
        os.path.join(VCF_DIR, '{indiv}.forimputation.vcf.gz')
    output:
        temp(os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.forimputation.recode.vcf.gz'))
    params:
        minimumdp = '{minDP}',
        prefix = os.path.join(VCF_DIR,'minDP{minDP}','{indiv}' + '.forimputation'),
        vcffile = os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.forimputation.recode.vcf'),
        logfile = os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.forimputation.log')
    shell:
        """
        {VCFTOOLS} --gzvcf {input} --min-meanDP {params.minimumdp} --recode --recode-INFO-all --out {params.prefix}
        rm {params.logfile}
        bgzip {params.vcffile}
        """


rule liftHg38ToHg19:
    input:
        os.path.join(VCF_DIR, 'minDP{minDP}', '{indiv}' + '.forimputation.recode.vcf.gz')
    params: 
        vcf_on_grch37 = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.vcf")
    output:
        middle=temp(os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.vcf.temp")),
        rejected = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.rejected.vcf"),
        vcf_gz = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.vcf.gz")
    shell:
        """
        zcat {input} |  grep "^#" > {output.middle}
        zcat {input} |  grep -v "^#" | awk "{{print "'"chr"'"\$0}}" >> {output.middle}
        {PICARD} LiftoverVcf I={output.middle} O={params.vcf_on_grch37} CHAIN={CHAIN_FN} REJECT={output.rejected} R={REF_FN}
        bgzip {params.vcf_on_grch37}
        """


rule format_vcf:
    input:
        vcf_gz = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.vcf.gz")
    params:
        os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.vcf.gz.tbi")
    output:
        temp(os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.header.txt"))
    shell:
        """
        rm -f {params}
        tabix -p vcf {input.vcf_gz}
        tabix -H {input.vcf_gz} > {output}
        """


rule split_chromosomes:
    input:
        vcf_gz = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.vcf.gz"),
        header = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.header.txt")
    params:
        chr_vcf = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}_chr{chr}.vcf")
    output:
        chr_vcf_gz = temp(os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}_chr{chr}.vcf.gz"))
    shell:
        """
        cat {input.header} > {params.chr_vcf}
        tabix {input.vcf_gz} chr{wildcards.chr} | sed 's/chr//g' >> {params.chr_vcf}
        bgzip {params.chr_vcf}
        tabix -p vcf {output.chr_vcf_gz}
        """

rule imputation:
    input:
        chr_vcf_gz = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}_chr{chr}.vcf.gz"),
        chr_REF_PANEL = "/work-zfs/abattle4/lab_data/imputation_reference_panel/{chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz"
    params:
        prefix = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}_chr{chr}.imputed"),
        info = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}_chr{chr}.imputed.info"),
        index_file = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}_chr{chr}.vcf.gz.tbi")
    output:
        temp(os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}_chr{chr}.imputed.dose.vcf.gz"))
    conda:
        "../envs/env_py37.yml"
    shell:
        """
        minimac4 --refHaps {input.chr_REF_PANEL} --haps {input.chr_vcf_gz} --prefix {params.prefix} --ChunkLengthMb 200 --ignoreDuplicates --format GT,DS,GP --minRatio 0.01
        rm {params.info}
        rm {params.index_file}
        """



rule merge_chrs:
    input:
        header = os.path.join(VCF_DIR, "minDP" + "{minDP}", "GRCh37", "{indiv}", "{indiv}.GRCh37.header.txt"),
        files = expand(os.path.join(VCF_DIR, "minDP" + "{{minDP}}", "GRCh37", "{{indiv}}", "{{indiv}}_chr{chr}.imputed.dose.vcf.gz"), chr = CHROM)
    output:
        vcf_file = os.path.join(BOWTIE_DIR, 'Imputation', 'minDP' + "{minDP}", "{indiv}", "{indiv}.imputed.GRCh37.vcf"),
    shell:
        """
        cat {input.header} | sed 's/chr//g'  > {output.vcf_file}
        zcat {input.files} | grep -v "#" >> {output.vcf_file}
        """


rule obtain_genotype_dosage:
    input:
        os.path.join(BOWTIE_DIR, 'Imputation', 'minDP' + "{minDP}", "{indiv}", "{indiv}.imputed.GRCh37.vcf")
    output:
        dosage = os.path.join(BOWTIE_DIR, 'Imputation', 'minDP' + "{minDP}", "{indiv}", "{indiv}.imputed.dosage.GRCh37.bed")
    conda:
        "../envs/env_py37.yml"
    shell:
        """
        {BCFTOOLS} query -f '%CHROM\t%POS\t[%GT\t%DS]\n' {input} | awk "{{print "'"chr"'"\$1,\$2,\$5 = \$2 + 1,\$3,\$4}}" > {output.dosage}
        """




LIFTOVER_DIR = "/work-zfs/abattle4/heyuan/tools/liftOver"

rule lift_dosage_to_GRCh38:
    input:
        dosage = os.path.join(BOWTIE_DIR, 'Imputation', 'minDP' + "{minDP}", "{indiv}", "{indiv}.imputed.dosage.GRCh37.bed"),
    output:
        lifted = os.path.join(BOWTIE_DIR, 'Imputation', 'minDP' + "{minDP}", "{indiv}", "{indiv}.imputed.dosage.GRCh38.bed"),
        unlifted = os.path.join(BOWTIE_DIR, 'Imputation', 'minDP' + "{minDP}", "{indiv}", "{indiv}.imputed.dosage.GRCh38.bed.unlifted")
    shell:
        """
        cd {LIFTOVER_DIR}
        ./liftOver {input.dosage} hg19ToHg38.over.chain.gz {output.lifted} {output.unlifted}
        """



