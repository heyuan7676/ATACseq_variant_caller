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
        expand(os.path.join(PEAK_DIR, '{indiv}' + '.id.bed'), indiv = INDIVS),
        os.path.join(PEAK_DIR, 'all.bed'),
        expand(os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed'), indiv = INDIVS),
        expand(os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed_matrix'), indiv = INDIVS)


'''Label your input BED files so that their IDs uniquely identify their intervals'''
rule label_bed_files:
    input:
        os.path.join(PEAK_DIR, '{indiv}' + '_peaks.narrowPeak')
    output:
        os.path.join(PEAK_DIR, '{indiv}' + '.id.bed')
    params:
        sample='{indiv}'
    shell:
        """        
        cut -f1-3 {input} | awk -vidx={{params.sample}} "{{ print \$0"\\t"idx; }}" > {output}
        """



'''Take the union of all these ID-tagged files with BEDOPS bedops'''
rule get_interval_incommon:
    input:
        expand(os.path.join(PEAK_DIR, '{indiv}' + '.id.bed'), indiv = INDIVS) 
    output:
        all = os.path.join(PEAK_DIR, 'all.bed'),
        union = os.path.join(PEAK_DIR, 'union-thresholded_level0.bed')
    shell:
        """
        bedops --everything {input} | bedmap --echo --echo-map-id-uniq --delim "\\t" - > {output.all}
        awk -vthreshold=3 "((split(\$5, ids, "'";"'")) >= threshold)" {output.all} | awk "{{print \$1,\$2,\$3,\$4}}" | sort -k1,1 -k2,2n | sed "s/ /	/g" > {output.union}
        """


'''Merge intervals that overlap more than 50%'''


UNION_PEAKS='union-peaks.bed'
rule peak_IDs:
    input:
        peaks = os.path.join(PEAK_DIR, 'union-thresholded_level8.bed')
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
        bedtools intersect -abam {input.bam_file} -b {input.union} -wo -bed | cut -d"	" -f13-16 | sort | uniq -c | awk "{{print \$5,\$1}}" | sort -k1,1n -k2,2n | sed "s/ /	/g" > {output}
        """

rule collect_peak_union:
    input:
        readscount = os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed'),
        peaks = os.path.join(PEAK_DIR, PEAK_UNION_fn)
    output:
        os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed_matrix')
    shell:
        """
        join -e0 -a 1 -a 2 -j 1 {input.peaks} -o auto {input.readscount} | awk "{{print \$5}}"> {output}
        """

rule obtain_peak_matrix:
    input:
        readscount = os.path.join(PEAK_DIR, '{indiv}.count.unionPeaks.bed_matrix'),
        peaks = os.path.join(PEAK_DIR, PEAK_UNION_fn)
    output:
        bychr = os.path.join(PEAK_DIR, "peak_by_sample_matrix_chr{chr}.txt")
    shell:
        """
        echo "PEAK CHR POS {INDIVS}" | tr "\\n" " " > {output.bychr}
        paste {input.peaks} {input.readscount} | awk "{{if((\$2 == "'"{wildcards.chr}"'")) print \$0}}" | sed 's/	/ /g' >> {output.bychr}
        """





