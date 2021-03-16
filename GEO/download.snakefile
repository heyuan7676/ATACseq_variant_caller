import os

'''Download dataset from GEO'''

rule download:
    output:
        reads1 = temp(os.path.join(FQ_DIR, '{indiv}' + '_1.fastq.gz')),
        reads2 = temp(os.path.join(FQ_DIR, '{indiv}' + '_2.fastq.gz')),
    params:
        sra_file = temp(os.path.join(FQ_DIR, '{indiv}' + '.sra')),
        outputdir = FQ_DIR,
        reads1 = temp(os.path.join(FQ_DIR, '{indiv}' + '_1.fastq')),
        reads2 = temp(os.path.join(FQ_DIR, '{indiv}' + '_2.fastq')),
    conda:
        "envs/env_py37.yml"
    shell:
        """
        rm -f {output}*
        {PREFETCH} {wildcards.indiv} -O {params.outputdir}
        cd {params.outputdir}
        {FADUMP} --split-files {params.sra_file}
        gzip {params.reads1}
        gzip {params.reads2}
        """

