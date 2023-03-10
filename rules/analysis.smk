

rule store_alignments_as_np_data:
    input:
        alignments="{path}/mapped.bam",
    output:
        "{path}/mapped.npz"
    shell:
        "numpy_alignments store -i {input.alignments} bam {output} -1"


rule get_accuracy_result:
    input:
        alignments=f"data/{parameters.until('n_threads')}/mapped.npz",
        truth=f"data/{parameters.until('n_reads')}/truth.npz"
    output:
        f"data/{parameters}/{{type, recall|one_minus_precision|f1_score}}.txt"
    params:
        allowed_bp_mismatch=50 #lambda wildcards: int(wildcards.read_length) // 5
    shell:
        "numpy_alignments get_correct_rates --report-type {wildcards.type} -m {wildcards.min_mapq} --allowed-bp-mismatch {params.allowed_bp_mismatch} {input.truth} {input.alignments} {wildcards.variant_filter} > {output}"


rule get_runtime:
    input:
        f"data/{parameters.until('n_threads')}/benchmark.csv"
    output:
        f"data/{parameters}/runtime.txt"
    shell:
        "cat {input} | tail -n 1 | cut -f 1 > {output}"



rule get_memory_usage:
    input:
        f"data/{parameters.until('n_threads')}/benchmark.csv"
    output:
        f"data/{parameters}/memory_usage.txt"
    shell:
        "cat {input} | tail -n 1 | cut -f 3 > {output}"


