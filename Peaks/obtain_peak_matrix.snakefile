import os

configfile: 'config.yaml'

BOWTIE_DIR = config['BOWTIE_DIR']
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks')
DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')

PEAK_UNION_fn = os.path.join(PEAK_DIR, 'union-peaks.bed')

INDIVS = glob_wildcards(os.path.join(PEAK_DIR, '{indiv}_peaks.narrowPeak'))
INDIVS = INDIVS[0]
CHROM = config['CHROM']
LEVELS = [0,1,2,3,4,5,6,7,8,9,10]

'''Collect peaks'''

rule all:
    input:
        expand(os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed'), indiv = INDIVS),
        expand(os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed_matrix'), indiv = INDIVS),
        expand(os.path.join(PEAK_DIR, "peak_by_sample_matrix_chr{chr}.txt"), chr = CHROM)

UNION_PEAKS='union-peaks.bed'
rule peak_IDs:
    input:
        peaks = os.path.join(PEAK_DIR, 'union-thresholded_level9.bed')
    output:
        union = os.path.join(PEAK_DIR, UNION_PEAKS)
    shell:
        """
        awk "{{print \$0,"'"Peak"'"NR}}" {input.peaks} | sed "s/ /	/g" > {output.union}
        """


# Count the number of reads overlapping with the peaks
rule count_reads:
    input:
        bam_file = os.path.join(DIR_FIRST_PASS, '{indiv}-clean.bam'),
        union = os.path.join(PEAK_DIR, UNION_PEAKS)
    output:
        os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed')
    shell:
        """
        bedtools intersect -abam {input.bam_file} -b {input.union} -wo -bed | cut -d"	" -f13-16 | sort | uniq -c | awk "{{print \$5,\$1}}" | sed "s/ /	/g" > {output}
        """

rule collect_peak_union:
    input:
        readscount = os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed'),
        peaks = os.path.join(PEAK_DIR, PEAK_UNION_fn)
    output:
        os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed_matrix')
    shell:
        """
        join -e0 -a 1 -a 2 -j 1 <(awk "{{print \$4}}" {input.peaks} | sort -k1,1) -o auto <(sort -k1,1 -k2,2n {input.readscount}) | awk "{{print \$2}}"> {output}
        """

rule obtain_peak_matrix:
    input:
        readscount = expand(os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed_matrix'), indiv = INDIVS),
        peaks = os.path.join(PEAK_DIR, PEAK_UNION_fn)
    output:
        bychr = os.path.join(PEAK_DIR, "peak_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "PEAK CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.peaks} {input.readscount} | awk "{{if((\$1 == "'"{wildcards.chr}"'")) print \$0}}" | sed 's/	/ /g' >> {output.bychr}
        """





