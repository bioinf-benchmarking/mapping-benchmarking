
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
        reads = "{data}/whole_genome_single_end/{config}/reads.fq.gz",
        idx = multiext("{data}/reference",".amb",".ann",".bwt",".pac",".sa"),
    output:
        "{data}/whole_genome_single_end/{config}/bwa/mapped.bam",
    params:
        extra=r"-R '@RG\tID:sample\tSM:sample'",
        sorting="none",# Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",# Can be 'queryname' or 'coordinate'.
        sort_extra="",# Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.21.2/bio/bwa/mem"


rule bwa_map_paired_end:
    input:
        reads = ["{data}/whole_genome_paired_end/{config}/reads" + n + ".fq.gz" for n in ("1", "2")],
        idx = multiext("{data}/reference",".amb",".ann",".bwt",".pac",".sa"),
    output:
        "{data}/whole_genome_paired_end/{config}/bwa/mapped.bam",
    params:
        extra=r"-R '@RG\tID:sample\tSM:sample'",
        sorting="none",# Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",# Can be 'queryname' or 'coordinate'.
        sort_extra="",# Extra args for samtools/picard.
    threads: 8
    wrapper:
        "v1.21.2/bio/bwa/mem"