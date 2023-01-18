

rule download_truth_vcf:
    output:
        "data/truth_variants/{individual}.vcf.gz"
    params:
        url=lambda wildcards: config["truth_datasets"][wildcards.individual]["vcf_url"]
    shell:
        "wget -O {output} {params.url}"



rule prepare_simulation:
    input:
        vcf="data/truth_variants/{individual}.vcf.gz",
        reference="data/reference_genomes/raw/{reference}.fa"

    output:
        coordinate_map="data/simulated_reads/{reference}/{individual}/coordinate_map_chromosome{chromosome}_haplotype{haplotype}.npz",
        haplotype_reference="data/simulated_reads/{reference}/{individual}/chromosome{chromosome}_haplotype{haplotype}_reference.fasta",
        #haplotype_reference="data/{dataset}/{truth_dataset}_chromosome{chromosome}_haplotype{haplotype}_reference.fasta",
        #haplotype_reference_fasta="data/{dataset}/{truth_dataset}_chromosome{chromosome}_haplotype{haplotype}_reference.fasta.fai",
    shell:
        "graph_read_simulator prepare_simulation --chromosome {wildcards.chromosome} --haplotype {wildcards.haplotype} "
        "--vcf {input.vcf} --reference {input.reference} -o data/simulated_reads/{wildcards.reference}/{wildcards.individual}/ "



"""
rule simulate_reads_for_chromosome_and_haplotype:
    input:
        coordinate_map="data/simulated_reads/{reference}/{individual}/coordinate_map_chromosome{chromosome}_haplotype{haplotype}.npz",
        haplotype_reference="data/simulated_reads/{reference}/{individual}/chromosome{chromosome}_haplotype{haplotype}_reference.fasta",
        haplotype_reference_fai="data/simulated_reads/{reference}/{individual}/chromosome{chromosome}_haplotype{haplotype}_reference.fasta.fai",
    output:
        "data/simulated_reads/{reference}/{individual}/raw_simulated_reads_chromosome{chromosome}_haplotype{haplotype}_coverage{coverage}.txt"
    shell:
        "graph_read_simulator simulate_reads -s 0.001 --deletion_prob 0.001 --insertion_prob 0.001 -D data/simulated_reads/{wildcards.reference}/{wildcards.individual}/ '{wildcards.chromosome} {wildcards.haplotype}' {wildcards.coverage} > {output}"
"""

rule simulate_reads_for_chromosome_and_haplotype:
    input:
        coordinate_map="data/simulated_reads/{reference}/{individual}/coordinate_map_chromosome{chromosome}_haplotype{haplotype}.npz",
        haplotype_reference="data/simulated_reads/{reference}/{individual}/chromosome{chromosome}_haplotype{haplotype}_reference.fasta",
        haplotype_reference_fai="data/simulated_reads/{reference}/{individual}/chromosome{chromosome}_haplotype{haplotype}_reference.fasta.fai",
    output:
        "data/simulated_reads/{reference}/{individual}/raw_simulated_reads_chromosome{chromosome}_haplotype{haplotype}_coverage{coverage}.txt"
    shell:
        "graph_read_simulator simulate_reads -s 0.001 --deletion_prob 0.001 --insertion_prob 0.001 -D data/simulated_reads/{wildcards.reference}/{wildcards.individual}/ '{wildcards.chromosome} {wildcards.haplotype}' {wildcards.coverage} > {output}"


def get_simulation_tmp_datasets(wildcards):
    haplotypes = [0, 1]
    chromosomes = config["simulations"][wildcards.dataset]["chromosomes"].split(",")
    genome_build = config["simulations"][wildcards.dataset]["genome"]
    individual = config["simulations"][wildcards.dataset]["individual"]
    coverage = str(config["simulations"][wildcards.dataset]["coverage"])

    files = []
    for chromosome in chromosomes:
        for haplotype in haplotypes:
            files.append("data/simulated_reads/" + genome_build + "/" + individual + "/raw_simulated_reads_chromosome" + chromosome + "_haplotype" + str(haplotype) + "_coverage" + coverage + ".txt")

    return files


rule simulate_reads:
    input:
        get_simulation_tmp_datasets
    output:
        reads="data/simulated_reads/{dataset}.fa",
        read_positions="data/simulated_reads/{dataset}.readpositions"
    shell:
        "cat {input} | graph_read_simulator assign_ids {output.read_positions} {output.reads}"

