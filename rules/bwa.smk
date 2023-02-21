
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
        idx = multiext(f"data/{parameters.until('dataset_size')}/reference", ".fa", ".fa.amb",".fa.ann",".fa.bwt",".fa.pac",".fa.sa"),
    output:
        f"data/{parameters.until('n_threads')(method='bwa')}/mapped.bam"
    benchmark:
        f"data/{parameters.until('n_threads')(method='bwa')}/benchmark.csv"
    threads: lambda wildcards: int(wildcards.n_threads)
    conda: "../envs/bwa.yml"
    shell:
        """
        bwa mem -t {wildcards.n_threads} -R "@RG\\tID:sample\\tSM:sample" {input.idx[0]} {input.reads} | samtools view -o {output} -
        """

