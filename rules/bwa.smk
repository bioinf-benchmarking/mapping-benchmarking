
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


rule bwa_map:
    input:
        reads = "data/simulated_reads/{dataset}/simulated_reads.fq.gz",
        idx = multiext("data/reference_genomes/{dataset}",".amb",".ann",".bwt",".pac",".sa"),
    output:
        "data/mapping/bwa/{dataset}.sam",
    params:
        extra=r"-R '@RG\tID:{dataset}\tSM:{dataset}'",
        sorting="none",# Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",# Can be 'queryname' or 'coordinate'.
        sort_extra="",# Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.21.2/bio/bwa/mem"