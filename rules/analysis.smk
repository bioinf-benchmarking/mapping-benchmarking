

rule store_alignments_as_np_data:
    input:
        alignments="{path}/{n_reads}/{method}/mapped.bam",
        #n_reads="data/simulated_reads/{dataset}.n_reads",
    output:
        "{path}/{n_reads,\d+}/{method}/mapped.npz"
    shell:
        "samtools view {input.alignments} | numpy_alignments store sam {output} {wildcards.n_reads}"


def get_result_files(wildcards):
    methods = config["runs"]["mapping"]["tools"]
    return [
        f"{wildcards.dataset}/{method}/mapped.npz" for method in methods
    ]


rule compare_read_mapping_against_truth:
    input:
        get_result_files,
        truth="{dataset}/truth.npz",
    output:
        "{dataset}/report.html"
    params:
        method_names=lambda wildcards: ','.join(config["runs"]["mapping"]["tools"]),
        input_files=lambda wildcards: ",".join(get_result_files(wildcards))
    shell:
        """
        numpy_alignments make_report -f {wildcards.dataset}/ --names="{params.method_names}" {input.truth} {params.input_files} purple,orange,blue
        """


rule make_mapping_report:
    output:
        "reports/mapping/report_{size}.html"



rule get_accuracy_result:
    input:
        alignments=f"data/{parameters.until('n_threads')}/mapped.npz",
        truth=f"data/{parameters.until('n_reads')}/truth.npz"
    output:
        f"data/{parameters}/{{type, recall|one_minus_precision|f1_score}}.txt"
    params:
        allowed_bp_mismatch=lambda wildcards: int(wildcards.read_length) // 5
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

