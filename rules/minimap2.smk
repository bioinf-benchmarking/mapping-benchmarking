

rule minimap2_map:
    input:
        query="data/simulated_reads/{dataset}.fa",
        target="data/reference_genomes/{dataset}.fa"
    output:
        sam="data/mapping/minimap2/{dataset}.sam",
    params:
        extra="-ax sr",# optional
        sorting="none",# optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",# optional: extra arguments for samtools/picard
    threads: 8
    wrapper:
        "v1.21.2/bio/minimap2/aligner"
