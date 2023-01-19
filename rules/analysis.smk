
"""
rule get_number_of_reads:
    input:
        "data/simulated_reads/{dataset}.readpositions"
    output:
        "data/simulated_reads/{dataset}.n_reads"
    shell:
        "wc -l {input} | cut -d ' ' -f 1 > {output}"
"""


"""
rule store_truth_positions_as_np_data:
    input:
        positions="data/simulated_reads/{dataset}.readpositions",
        n_reads="data/simulated_reads/{dataset}.n_reads",
    output:
        positions="data/simulated_reads/{dataset}.truth.npz",
    params:
        n_reads=lambda wildcards: open("data/simulated_reads/" + wildcards.dataset + ".n_reads").read().strip()
    shell:
        "cat {input.positions} | numpy_alignments store truth {output} {params.n_reads}"
"""


rule store_alignments_as_np_data:
    input:
        alignments="data/mapping/{method}/{dataset}.sam",
        #n_reads="data/simulated_reads/{dataset}.n_reads",
    output:
        "data/mapping/{method}/{dataset}.npz"
    params:
        n_reads = lambda wildcards: config["simulations"][wildcards.dataset]["n_reads"]
        #n_reads=lambda wildcards: open("data/simulated_reads/" + wildcards.dataset + ".n_reads").read().strip()
    shell:
        "cat {input.alignments} | numpy_alignments store sam {output} {params.n_reads}"


def get_result_files(wildcards):
    methods = config["runs"]["mapping"]["tools"]
    return [
        f"data/mapping/{method}/{wildcards.dataset}.npz" for method in methods
    ]


rule compare_read_mapping_against_truth:
    input:
        get_result_files,
        truth="data/simulated_reads/{dataset}/truth.npz",
        #bwa="data/mapping/bwa/{dataset}.npz",
        #strobealign="data/mapping/strobealign/{dataset}.npz",
        #minimap2="data/mapping/minimap2/{dataset}.npz",
    output:
        "reports/mapping/{dataset}/report.html"
    params:
        method_names=lambda wildcards: ','.join(config["runs"]["mapping"]["tools"]),
        input_files=lambda wildcards: ",".join(get_result_files(wildcards))
    shell:
        """
        numpy_alignments make_report -f reports/mapping/{wildcards.dataset}/ --names="{params.method_names}" {input.truth} {params.input_files} purple,orange,blue
        """


rule make_mapping_report:
    output:
        "reports/mapping/report_{size}.html"




