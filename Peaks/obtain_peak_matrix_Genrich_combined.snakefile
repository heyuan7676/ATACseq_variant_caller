import os

configfile: 'config.yaml'

BOWTIE_DIR = config['BOWTIE_DIR']
PEAK_DIR = os.path.join(BOWTIE_DIR, 'Peaks_Genrich')
Peak_OUTPUT_DIR = os.path.join(PEAK_DIR, 'combined')

DIR_FIRST_PASS = os.path.join(BOWTIE_DIR, 'first_pass_bqsr')
if not os.path.isdir(Peak_OUTPUT_DIR):
    os.makedirs(Peak_OUTPUT_DIR)

PEAK_UNION_fn = os.path.join(Peak_OUTPUT_DIR, 'union-peaks_combined.bed')

INDIVS = glob_wildcards(os.path.join(PEAK_DIR, '{indiv}_peaks.Genrich.narrowPeak'))
INDIVS = INDIVS[0]
INDIVS = [x for x in INDIVS if x.startswith('HG')]
INDIVS.sort()

CHROM = config['CHROM']

'''Collect peaks'''

rule all:
    input:
        expand(os.path.join(Peak_OUTPUT_DIR, '{indiv}.count.unionPeaks.bed'), indiv = INDIVS),
        expand(os.path.join(Peak_OUTPUT_DIR, '{indiv}.count.unionPeaks.bed_matrix'), indiv = INDIVS),
        expand(os.path.join(Peak_OUTPUT_DIR, "peak_by_sample_matrix_chr{chr}.txt"), chr = CHROM)

rule peak_IDs:
    input:
        peaks = os.path.join(PEAK_DIR, 'combined_peaks.Genrich.narrowPeak')
    output:
        union = PEAK_UNION_fn
    shell:
        """
        awk "{{print \$0,"'"Peak"'"NR}}" {input.peaks} | sort -k4,4 | sed "s/ /	/g  " > {output.union}
        """


# Count the number of reads overlapping with the peaks
rule count_reads:
    input:
        bam_file = os.path.join(DIR_FIRST_PASS, '{indiv}-clean.bam'),
        union = PEAK_UNION_fn
    output:
        temp(os.path.join(Peak_OUTPUT_DIR, '{indiv}.count.unionPeaks.bed'))
    shell:
        """
        bedtools intersect -abam {input.bam_file} -b {input.union} -wo -bed |  sed "s/\/1//g" | sed "s/\/2//g" | awk "{{print \$4,\$16}}" | sort | uniq | awk "{{print \$2}}" | sort | uniq -c | awk "{{print \$2,\$1}}" | sort -k1,1 -k2,2n | sed "s/ /	/g " > {output}
        """

rule collect_peak_union:
    input:
        readscount = os.path.join(Peak_OUTPUT_DIR, '{indiv}.count.unionPeaks.bed'),
        peaks = os.path.join(Peak_OUTPUT_DIR, PEAK_UNION_fn)
    output:
        temp(os.path.join(Peak_OUTPUT_DIR, '{indiv}.count.unionPeaks.bed_matrix'))
    shell:
        """
        join -e0 -a 1 -a 2 -j 1 <(awk "{{print \$4}}" {input.peaks}) -o auto <(cat {input.readscount}) | awk "{{print \$2}}"> {output}
        """

rule obtain_peak_matrix:
    input:
        readscount = expand(os.path.join(Peak_OUTPUT_DIR, '{indiv}.count.unionPeaks.bed_matrix'), indiv = INDIVS),
        peaks = PEAK_UNION_fn
    output:
        bychr = os.path.join(Peak_OUTPUT_DIR, "peak_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "CHR START END PEAK {INDIVS}" | tr "\\n" " " > {output.bychr}
        echo "" >> {output.bychr}
        paste {input.peaks} {input.readscount} | awk "{{if((\$1 == "'"{wildcards.chr}"'")) print \$0}}" | sed 's/	/ /g' >> {output.bychr}
        """


