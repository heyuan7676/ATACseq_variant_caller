import os

configfile: 'config.yaml'

PICARD = config['PICARD']
BCFTOOLS = config['BCFTOOLS']
VCFTOOLS = config['VCFTOOLS']
minDP = config['minDP']

DATA_DIR = config['DATA_DIR']
VCF_DIR = os.path.join(DATA_DIR, 'VCF_files')

GRCH37_DIR = os.path.join(VCF_DIR, 'GRCh37')
if not os.path.isdir(GRCH37_DIR):
    os.makedirs(GRCH37_DIR)

IMPUTATION_DIR = os.path.join(DATA_DIR, "Imputation", "minDP" + str(minDP))
if not os.path.isdir(IMPUTATION_DIR):
    os.makedirs(IMPUTATION_DIR)

INDIVS = glob_wildcards(os.path.join(VCF_DIR, '{indiv}.filtered.recode.vcf.gz'))
INDIVS = INDIVS[0]
INDIVS = INDIVS[:1]

for s in INDIVS:
    S_DIR = os.path.join(GRCH37_DIR, s)
    if not os.path.isdir(S_DIR):
        os.makedirs(S_DIR)
    S_DIR = os.path.join(IMPUTATION_DIR, s)
    if not os.path.isdir(S_DIR):
        os.makedirs(S_DIR)

CHROM = config['CHROM']
CHROM = ['22']

'''Collect genoyptes'''

wildcard_constraints:
    chr="\d+"

rule all:
    input:
        #expand(os.path.join(GRCH37_DIR, "{indiv}.filtered.recode.GRCh37.vcf.gz"), indiv = INDIVS),
        #expand(os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.withINFO.vcf.gz"), indiv = INDIVS, chr = CHROM)
        expand(os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.imputed.GRCh37.minDP" + minDP + ".dose.vcf.gz"), indiv = INDIVS, chr = CHROM),
        expand(os.path.join(IMPUTATION_DIR, "{indiv}", "chr{chr}.imputed.GRCh38.genotype.txt"), indiv = INDIVS, chr = CHROM)


CHAIN_FN = "/work-zfs/abattle4/heyuan/tools/liftOver/hg38ToHg19.over.chain.gz"
REF_FN = "/work-zfs/abattle4/heyuan/database/GRCh37_reference/hg19.fa"

rule liftHg38ToHg19:
    input:
        vcf_file = os.path.join(VCF_DIR, "{indiv}.filtered.recode.vcf.gz")
    params: 
        vcf_file = os.path.join(VCF_DIR, "{indiv}.filtered.recode.vcf"),
        vcf_on_grch37 = os.path.join(GRCH37_DIR, "{indiv}.filtered.recode.GRCh37.vcf")
    output:
        middle=temp(os.path.join(VCF_DIR, "{indiv}.filtered.recode.vcf.temp")),
        rejected = os.path.join(GRCH37_DIR, "{indiv}.filtered.recode.GRCh37.rejected.vcf"),
        vcf_gz = os.path.join(GRCH37_DIR, "{indiv}.filtered.recode.GRCh37.vcf.gz")
    shell:
        """
        gunzip {input.vcf_file}
        cat {params.vcf_file} |  grep "^#" > {output.middle}
        cat {params.vcf_file} |  grep -v "^#" | awk "{{print "'"chr"'"\$0}}" >> {output.middle}
        {PICARD} LiftoverVcf I={output.middle} O={params.vcf_on_grch37} CHAIN={CHAIN_FN} REJECT={output.rejected} R={REF_FN}
        bgzip {params.vcf_file}
        bgzip {params.vcf_on_grch37}
        """


rule format_vcf:
    input:
        vcf_gz = os.path.join(GRCH37_DIR, "{indiv}.filtered.recode.GRCh37.vcf.gz")
    output:
        os.path.join(GRCH37_DIR, "{indiv}", "{indiv}.filtered.recode.GRCh37.header.txt")
    shell:
        """
        tabix -p vcf {input.vcf_gz}
        tabix -H {input.vcf_gz} > {output}
        """


rule split_chromosomes:
    input:
        vcf_gz = os.path.join(GRCH37_DIR, "{indiv}.filtered.recode.GRCh37.vcf.gz"),
        header = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}.filtered.recode.GRCh37.header.txt")
    params:
        chr_vcf = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.withINFO.vcf")
    output:
        chr_vcf_gz = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.withINFO.vcf.gz"),
    shell:
        """
        cat {input.header} > {params.chr_vcf}
        tabix {input.vcf_gz} chr{wildcards.chr} | sed 's/chr//g' | sed '/	<NON_REF>/d' | sed 's/,<NON_REF>//g' >> {params.chr_vcf}
        bgzip {params.chr_vcf}
        tabix -p vcf {output.chr_vcf_gz}
        """


rule filterVCF_minDP:
    input:
        os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.withINFO.vcf.gz")
    params:
        minimumdp = minDP,
        prefix = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.minDP"+ minDP),
        chr_vcf = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.minDP"+ minDP + ".recode.vcf"),
        tempfn = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.minDP"+ minDP + ".recode.vcf.temp")
    output:
        chr_vcf_gz = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.minDP" + minDP + ".recode.vcf.gz"),
    shell:
        """
        {VCFTOOLS} --gzvcf {input} --min-meanDP {params.minimumdp} --recode --recode-INFO-all --out {params.prefix}
        bgzip {params.chr_vcf}
        tabix -f -p vcf {output.chr_vcf_gz}
        {BCFTOOLS} annotate -x INFO {output.chr_vcf_gz} |  {BCFTOOLS} annotate -x FORMAT > {params.tempfn}
        rm {output.chr_vcf_gz}
        mv {params.tempfn} {params.chr_vcf}
        bgzip {params.chr_vcf}
        tabix -f -p vcf {output.chr_vcf_gz}
        """


rule imputation:
    input:
        chr_vcf_gz = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.minDP" + minDP + ".recode.vcf.gz"),
        chr_ref_panel = "/work-zfs/abattle4/lab_data/imputation_reference_panel/{chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz"
    params:
        prefix = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.imputed.GRCh37.minDP" + minDP)
    output:
        os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.imputed.GRCh37.minDP" + minDP + ".dose.vcf.gz")
    conda:
        "../envs/env_py37.yml"
    shell:
        'minimac4 --refHaps {input.chr_ref_panel} --haps {input.chr_vcf_gz} --prefix {params.prefix} --ChunkLengthMb 100 --ignoreDuplicates --format GT,DS,GP'





rule obtain_genotype:
    input:
        os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.imputed.GRCh37.minDP" + minDP+ ".dose.vcf.gz")
    params:
        intermediate = os.path.join(GRCH37_DIR, "{indiv}", "{indiv}_chr{chr}.imputed.GRCh37.minDP" + minDP + ".dose.vcf")
    output:
        os.path.join(IMPUTATION_DIR, "{indiv}", "chr{chr}.imputed.GRCh37.genotype.txt")
    conda:
        "../envs/env_py37.yml"
    shell:
        """
        gunzip -f {input}
        vcf-to-tab < {params.intermediate} > {output}
        """




LIFTOVER_DIR = "/work-zfs/abattle4/heyuan/tools/liftOver"
rule lift_to_GRCh38:
    input:
        os.path.join(IMPUTATION_DIR, "{indiv}", "chr{chr}.imputed.GRCh37.genotype.txt")
    output:
        bedfile = temp(os.path.join(IMPUTATION_DIR, "{indiv}", "chr{chr}.imputed.GRCh37.genotype.bed")),
        lifted = os.path.join(IMPUTATION_DIR, "{indiv}", "chr{chr}.imputed.GRCh38.genotype.txt"),
    shell:
        """
        awk "{{print "'"chr"'"\$1,\$2,\$3 = \$2 + 1,\$4}}" {input} | sed "1d" > {output.bedfile}
        cd {LIFTOVER_DIR}
        ./liftOver {output.bedfile} hg19ToHg38.over.chain.gz {output.lifted} unlifted.bed
        """ 

