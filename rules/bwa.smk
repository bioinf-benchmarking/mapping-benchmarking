from mapping_benchmarking.parameter_config import MappedReads, ReferenceGenome, Reads, SingleEndReads

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
        #reads = get_input_reads,
        reads = SingleEndReads.as_input(),
        idx = ReferenceGenome.path(file_ending=[".fa", ".fa.amb",".fa.ann",".fa.bwt",".fa.pac",".fa.sa"])
        #idx = multiext(f"data/{reference_genome}/reference", ".fa", ".fa.amb",".fa.ann",".fa.bwt",".fa.pac",".fa.sa")
    output:
        reads=MappedReads.as_output(method='bwa')
        #f"data/{reference_genome}/{{config}}/bwa/{{n_threads}}/mapped.bam"
    benchmark:
        #f"data/{parameters.until('n_threads')(method='bwa')}/benchmark.csv"
        #f"data/{parameters.until('dataset_size')}/{{config}}/bwa/{{n_threads}}/benchmark.csv"
        MappedReads.as_output(method='bwa', file_ending=".benchmark.csv")
        #f"data/{reference_genome}/{{config}}/bwa/{{n_threads}}/benchmark.csv"
    threads: lambda wildcards: int(wildcards.n_threads)
    conda: "../envs/bwa.yml"
    shell:
        """
        bwa mem -t {wildcards.n_threads} -R "@RG\\tID:sample\\tSM:sample" {input.idx[0]} {input.reads} | samtools view -o {output} -
        """

