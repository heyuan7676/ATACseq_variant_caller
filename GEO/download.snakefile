import os

'''Download dataset from GEO'''

rule download:
    output:
        temp(os.path.join(FQ_DIR, '{indiv}' + '.sra'))
    params:
        outputdir = FQ_DIR
    conda:
        "envs/env_py37.yml"
    shell:
        """
        rm -f {output}*
        {PREFETCH} {wildcards.indiv} -O {params.outputdir}
        """


rule transform_to_fa:
    input:
        os.path.join(FQ_DIR, '{indiv}' + '.sra')
    output:
        reads1 = temp(os.path.join(FQ_DIR, '{indiv}' + '_1.fastq')),
        reads2 = temp(os.path.join(FQ_DIR, '{indiv}' + '_2.fastq')),
    params:
        outputdir = FQ_DIR
    shell:
        """
        cd {params.outputdir}
        {FADUMP} --split-files {input}
        """
