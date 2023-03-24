from mapping_benchmarking.config import *


rule bwa_index:
    input:
        "{genome}.fa",
    output:
        idx=multiext("{genome}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.2/bio/bwa/index"


rule bwa_map:
    input:
        reads = get_input_reads,
        idx = ReferenceGenome.path(file_ending=[".fa", ".fa.amb",".fa.ann",".fa.bwt",".fa.pac",".fa.sa"])
    output:
        reads=GenericMappedReads.as_output(method='bwa')
    benchmark:
        GenericMappedReads.as_output(method='bwa', file_ending=".benchmark.csv")
    threads: lambda wildcards: int(wildcards.n_threads)
    conda: "../envs/bwa.yml"
    shell:
        """
        bwa mem -t {wildcards.n_threads} -R "@RG\\tID:sample\\tSM:sample" {input.idx[0]} {input.reads} | samtools view -o {output} -
        """

