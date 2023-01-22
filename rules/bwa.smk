
rule bwa_index:
    input:
        "{genome}.fa",
    output:
        idx=multiext("{genome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/{genome}.log",
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.2/bio/bwa/index"


rule bwa_map_single_end:
    input:
        reads = "data/simulated_reads/{individual}/{config}/reads.fq.gz",
        idx = multiext("data/simulated_reads/{individual}/reference_genome",".amb",".ann",".bwt",".pac",".sa"),
    output:
        "data/simulated_reads/{individual,[a-z0-9_]+}/{config}/bwa/mapped.bam",
    params:
        extra=r"-R '@RG\tID:{individual}\tSM:{individual}'",
        sorting="none",# Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",# Can be 'queryname' or 'coordinate'.
        sort_extra="",# Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.21.2/bio/bwa/mem"