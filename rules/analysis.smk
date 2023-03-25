from mapping_benchmarking.config import WholeGenomeMappedReads, MappingAccuracy, GenericReads, GenericMappedReads


rule store_alignments_as_np_data:
    input:
        alignments=GenericMappedReads.path(file_ending=".bam"),
    output:
        GenericMappedReads.path(file_ending=".npz")
    shell:
        "numpy_alignments store -i {input.alignments} bam {output} -1"


rule get_accuracy_result:
    input:
        alignments=WholeGenomeMappedReads.path(file_ending=".npz"),
        truth=WholeGenomeReads.path(file_ending="/truth.npz")
    output:
        MappingAccuracy.path()
    params:
        allowed_bp_mismatch=50
    shell:
        "numpy_alignments get_correct_rates --report-type {wildcards.accuracy_type} -m {wildcards.min_mapq} --allowed-bp-mismatch {params.allowed_bp_mismatch} {input.truth} {input.alignments} {wildcards.variant_filter} > {output}"


rule get_runtime:
    input:
        WholeGenomeMappedReads.as_output(file_ending=".benchmark.csv")
        #f"data/{parameters.until('n_threads')}/benchmark.csv"
    output:
        Runtime.path()
        #f"data/{parameters}/runtime.txt"
    shell:
        "cat {input} | tail -n 1 | cut -f 1 > {output}"



rule get_memory_usage:
    input:
        WholeGenomeMappedReads.as_output(file_ending=".benchmark.csv")
        #f"data/{parameters.until('n_threads')}/benchmark.csv"
    output:
        MemoryUsage.path()
        #f"data/{parameters}/memory_usage.txt"
    shell:
        "cat {input} | tail -n 1 | cut -f 3 > {output}"


