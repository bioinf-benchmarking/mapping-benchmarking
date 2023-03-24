from mapping_benchmarking.config import *


rule bowtie2_index:
    input:
        ref="{genome}.fa",
    output:
        multiext(
            "{genome}.fa",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    params:
        extra="--large-index",  # optional parameters
    threads: 4
    wrapper:
        "v1.23.1/bio/bowtie2/build"



rule bowtie2_map:
    input:
        reads = get_input_reads,
        idx = ReferenceGenome.path(file_ending=[
            ".fa",
            ".fa.1.bt2l",
            ".fa.2.bt2l",
            ".fa.3.bt2l",
            ".fa.4.bt2l",
            ".fa.rev.1.bt2l",
            ".fa.rev.2.bt2l"
            ]
        )
    output:
        reads=GenericMappedReads.as_output(method="bowtie2")
        #f"data/{parameters.until('n_threads')(method='bowtie2', read_type='whole_genome_single_end')}/mapped.bam"
        #f"data/{reference_genome}/{{config}}/bowtie2/{{n_threads}}/mapped.bam"
    benchmark:
        #f"data/{parameters.until('n_threads')(method='bowtie2', read_type='whole_genome_single_end')}/benchmark.csv"
        GenericMappedReads.as_output(method='bowtie2', file_ending=".benchmark.csv")
        #f"data/{reference_genome}/{{config}}/bowtie2/{{n_threads}}/benchmark.csv"
    threads: lambda wildcards: int(wildcards.n_threads)
    conda: "../envs/bowtie2.yml"
    params:
        input_reads_string=lambda wildcards, input, output: f"-U {input.reads}" if isinstance(input.reads, str) else f"-1 {input.reads[0]} -2 {input.reads[2]}"
    shell:
        """
        bowtie2 --rg-id sample --rg "SM:sample" -x {input.idx[0]} -p {wildcards.n_threads} {params.input_reads_string} | samtools view -b -h - > {output} 
        """

