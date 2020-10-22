import os

rule sub_sample:
    input:
        #os.path.join(BOWTIE_DIR, '{indiv}' + SUFFIX + '.bam')
        os.path.join(DIR_FIRST_PASS, '{indiv}' + '-clean.bam')
    output:
        os.path.join(SUBSAMPLE_DIR, '{indiv}' + '-clean_subsampled.bam')
    params:
        prop=Proportion
    shell:
        '{SAMTOOLS} view -s 10.{params.prop} -b {input} > {output}'
