

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
        alignments="{config}/mapped.npz",
        truth=lambda wildcards: "/".join(wildcards.config.split("/")[:-2])  + "/truth.npz"  # hacky, get truth file two dirs down
    output:
        "{config}/{mapq}/{variant_filter}/{type}.txt"  # type is recall, one_minus_precision or f1_score
    shell:
        "numpy_alignments get_correct_rates --report-type {wildcards.type} -m {wildcards.mapq} {input.truth} {input.alignments} {wildcards.variant_filter} > {output}"



rule get_runtime:
    input:
        "{path}/benchmark.csv"
    output:
        "{path}/runtime.txt"
    shell:
        "cat {input} | tail -n 1 | cut -f 1 > {output}"