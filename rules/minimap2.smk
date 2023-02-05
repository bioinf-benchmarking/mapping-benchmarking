



rule minimap_single_end:
    input:
        query = "{data}/whole_genome_single_end/{config}/reads.fq.gz",
        target = "{data}/reference.fa"
    output:
        "{data}/whole_genome_single_end/{config}/minimap/{n_threads}/mapped.bam",
    params:
        extra = "-ax sr",  # optional
        sorting = "none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra = "",  # optional: extra arguments for samtools/picard

    benchmark:
        "{data}/whole_genome_single_end/{config}/minimap/{n_threads}/benchmark.csv",
    threads: lambda wildcards: int(wildcards.n_threads)
    wrapper:
        "v1.21.2/bio/minimap2/aligner"


rule minimap_paired_end:
    input:
        query = ["{data}/whole_genome_paired_end/{config}/reads" + n + ".fq.gz" for n in ("1", "2")],
        target = "{data}/reference.fa"
    output:
        "{data}/whole_genome_paired_end/{config}/minimap/{n_threads}/mapped.bam",
    params:
        extra="-ax sr",# optional
        sorting="none",# optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",# optional: extra arguments for samtools/picard
    threads: lambda wildcards: int(wildcards.n_threads)
    wrapper:
        "v1.21.2/bio/minimap2/aligner"
