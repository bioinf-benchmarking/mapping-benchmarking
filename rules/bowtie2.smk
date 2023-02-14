
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
    threads: 8
    wrapper:
        "v1.23.1/bio/bowtie2/build"



rule bowtie2_map_single_end:
    input:
        reads = single_end_reads,
        idx=multiext(
            f"data/{parameters.until('dataset_size')}/reference",
            ".fa",
            ".fa.1.bt2l",
            ".fa.2.bt2l",
            ".fa.3.bt2l",
            ".fa.4.bt2l",
            ".fa.rev.1.bt2l",
            ".fa.rev.2.bt2l",
        ),
    output:
        f"data/{parameters.until('n_threads')(method='bowtie2', read_type='whole_genome_single_end')}/mapped.bam"
    benchmark:
        f"data/{parameters.until('n_threads')(method='bowtie2', read_type='whole_genome_single_end')}/benchmark.csv"
    threads: lambda wildcards: int(wildcards.n_threads)
    conda: "../envs/bowtie2.yml"
    shell:
        """
        bowtie2 --rg-id sample --rg "SM:sample" -x {input.idx[0]} -p {wildcards.n_threads} -U {input.reads} | samtools view -b -h - > {output} 
        """


rule bowtie2_map_paired_end:
    input:
        reads=paired_end_reads,
        idx=multiext(
            f"data/{parameters.until('dataset_size')}/reference",
            ".fa",
            ".fa.1.bt2l",
            ".fa.2.bt2l",
            ".fa.3.bt2l",
            ".fa.4.bt2l",
            ".fa.rev.1.bt2l",
            ".fa.rev.2.bt2l",
        ),
    output:
        f"data/{parameters.until('n_threads')(method='bowtie2', read_type='whole_genome_paired_end')}/mapped.bam"
    benchmark:
        f"data/{parameters.until('n_threads')(method='bowtie2', read_type='whole_genome_paired_end')}/benchmark.csv"
    threads: lambda wildcards: int(wildcards.n_threads)
    conda: "../envs/bowtie2.yml"
    shell:
        """
        bowtie2 --rg-id sample --rg "SM:sample" -x {input.idx[0]} -p {wildcards.n_threads} -1 {input.reads[0]} -2 {input.reads[1]} | samtools view -b -h - > {output} 
        """


