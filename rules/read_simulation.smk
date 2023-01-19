

rule download_truth_vcf:
    output:
        "data/truth_variants/{individual}.vcf.gz"
    params:
        url=lambda wildcards: config["truth_datasets"][wildcards.individual]["vcf_url"]
    shell:
        "wget -O {output} {params.url}"



rule make_chromosome_haplotype_sequence_for_simulation:
    input:
        vcf="data/truth_variants/{individual}.vcf.gz",
        reference="data/reference_genomes/raw/{reference}.fa"

    output:
        coordinate_map="data/simulated_reads/{reference}/{individual}/coordinate_map_chromosome{chromosome}_haplotype{haplotype}.npz",
        haplotype_reference="data/simulated_reads/{reference}/{individual}/chromosome{chromosome}_haplotype{haplotype}_reference.fasta",
    shell:
        "graph_read_simulator prepare_simulation --chromosome {wildcards.chromosome} --haplotype {wildcards.haplotype} "
        "--vcf {input.vcf} --reference {input.reference} -o data/simulated_reads/{wildcards.reference}/{wildcards.individual}/ "


def get_input_files_for_haplotype_sequence(wildcards):

    chromosomes = config["simulations"][wildcards.dataset]["chromosomes"].split(",")
    genome_build = config["simulations"][wildcards.dataset]["genome"]
    individual = config["simulations"][wildcards.dataset]["individual"]

    return [
        f"data/simulated_reads/{genome_build}/{individual}/chromosome{chromosome}_haplotype{wildcards.haplotype}_reference.fasta"
        for chromosome in chromosomes
    ]

rule make_haplotype_sequence_for_simulation:
    input:
        get_input_files_for_haplotype_sequence
    output:
        "data/simulated_reads/{dataset}/haplotype{haplotype}.fa",
    shell:
        """
        cat {input} > {output}
        """

def get_haplotype_coordinate_maps(wildcards):

    chromosomes = config["simulations"][wildcards.dataset]["chromosomes"].split(",")
    genome_build = config["simulations"][wildcards.dataset]["genome"]
    individual = config["simulations"][wildcards.dataset]["individual"]

    return [
        f"data/simulated_reads/{genome_build}/{individual}/coordinate_map_chromosome{chromosome}_haplotype{wildcards.haplotype}.npz"
        for chromosome in chromosomes
    ]

rule merge_haplotype_coordinate_maps:
    input:
        get_haplotype_coordinate_maps
    output:
        "data/simulated_reads/{dataset}/coordinate_map_haplotype{haplotype}.npz",
    run:
        from shared_memory_wrapper import to_file
        from graph_read_simulator.simulation import CoordinateMap, MultiChromosomeCoordinateMap
        chromosomes = config["simulations"][wildcards.dataset]["chromosomes"].split(",")
        data = {}
        for chromosome, file in zip(chromosomes, input):
            data[chromosome] = CoordinateMap.from_file(file)
        to_file(MultiChromosomeCoordinateMap(data), output[0])

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
        haplotype_reference="data/simulated_reads/{dataset}/haplotype{haplotype}.fa"
    output:
        simulated_reads="data/simulated_reads/{dataset}/simulated_reads_haplotype{haplotype}.fq.gz",
        truth_positions="data/simulated_reads/{dataset}/simulated_reads_haplotype{haplotype,\d+}.sam"
    conda:
        "../envs/mason.yml"
    params:
        n_reads=lambda wildcards: config["simulations"][wildcards.dataset]["n_reads"] // 2  # divide by 2 for haplotype
    shell:
        """
        mason_simulator -ir {input.haplotype_reference} -n {params.n_reads} -o {output.simulated_reads} -oa {output.truth_positions}
        """


rule change_truth_alignments_to_reference_coordinates:
    input:
        truth_positions = "data/simulated_reads/{dataset}/simulated_reads_haplotype{haplotype}.sam",
        coordinate_map = "data/simulated_reads/{dataset}/coordinate_map_haplotype{haplotype}.npz",
    output:
        "data/simulated_reads/{dataset}/simulated_reads_haplotype{haplotype,\d+}.reference_coordinates.sam"
    run:
        from shared_memory_wrapper import from_file
        coordinate_map = from_file(input.coordinate_map)
        out_file = open(output[0], "w")
        for line in open(input.truth_positions):
            if line.startswith("@"):
                out_file.write(line)
            else:
                l = line.split("\t")
                chromosome = l[2]
                position = int(l[3])
                ref_position = coordinate_map.convert(chromosome, position)
                l[3] = str(ref_position)
                out_file.write("\t".join(l))


rule merge_simulated_reads_for_haplotypes:
    input:
        haplotype0="data/simulated_reads/{dataset}/simulated_reads_haplotype0.fq.gz",
        haplotype1="data/simulated_reads/{dataset}/simulated_reads_haplotype1.fq.gz",
    output:
        "data/simulated_reads/{dataset}/simulated_reads.fq.gz"
    shell:
        "zcat {input} | python scripts/assign_ids_to_fq.py | gzip > {output} "


rule merge_truth_alignments_for_haplotypes:
    input:
        haplotype0 = "data/simulated_reads/{dataset}/simulated_reads_haplotype0.reference_coordinates.sam",
        haplotype1 = "data/simulated_reads/{dataset}/simulated_reads_haplotype1.reference_coordinates.sam",
    output:
        "data/simulated_reads/{dataset}/simulated_reads.sam"
    shell:
        "cat {input} | python scripts/assign_ids_to_sam.py  > {output} "


rule store_truth_alignments:
    input:
        "data/simulated_reads/{dataset}/simulated_reads.sam"
    output:
        "data/simulated_reads/{dataset}/truth.npz"
    params:
        n_reads = lambda wildcards: config["simulations"][wildcards.dataset]["n_reads"]
    shell:
        """
        cat {input} | numpy_alignments store sam {output} {params.n_reads}
        """


"""
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

"""