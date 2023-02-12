
rule bowtie2_index:
    input:
        ref="{genome}.fa",
    output:
        multiext(
            "{genome}",
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
        sample = "{data}/whole_genome_single_end/{config}/reads.fq.gz",
        #sample=["reads/{sample}.1.fastq", "reads/{sample}.2.fastq"],
        idx=multiext(
            "{data}/reference",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    output:
        "{data}/whole_genome_single_end/{config}/bowtie2/{n_threads}/mapped.bam",
    benchmark:
        "{data}/whole_genome_single_end/{config}/bowtie2/{n_threads}/benchmark.csv",
    params:
        extra="--rg-id sample --rg SM:sample",  # optional parameters
    threads: 2  # Use at least two threads
    wrapper:
        "v1.23.1/bio/bowtie2/align"


rule bowtie2_map_paired_end:
    input:
        sample = ["{data}/whole_genome_paired_end/{config}/reads" + n + ".fq.gz" for n in ("1", "2")],
        idx=multiext(
            "{data}/reference",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    output:
        "{data}/whole_genome_paired_end/{config}/bowtie2/{n_threads}/mapped.bam",
    params:
        extra="--rg-id sample --rg SM:sample",  # optional parameters
    benchmark:
        "{data}/whole_genome_paired_end/{config}/bowtie2/{n_threads}/benchmark.csv",
    threads: 2  # Use at least two threads
    wrapper:
        "v1.23.1/bio/bowtie2/align"
