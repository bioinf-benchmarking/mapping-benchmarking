

rule download_truth_vcf:
    output:
        "data/{genome_build}/{individual}/variants.vcf.gz"
    params:
        url=lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.individual]["vcf_url"]
    shell:
        "wget -O {output} {params.url}"


rule download_truth_regions:
    output:
        "data/{genome_build}/{individual}/truth_regions.bed"
    params:
        url=lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.individual]["truth_regions"]
    conda:
        "../envs/bcftools.yml"
    shell:
        "wget -O {output} {params.url}"


def get_individual_properties(wildcards, property):
    return config["individuals"][wildcards.dataset][property]


rule make_chromosome_haplotype_sequence_for_simulation:
    input:
        vcf="data/{genome_build}/{individual}/variants.vcf.gz",
        reference="data/{genome_build}/reference.fa"

    output:
        coordinate_map="data/{genome_build}/{individual}/coordinate_map_chromosome{chromosome}_haplotype{haplotype}.npz",
        haplotype_reference="data/{genome_build}/{individual}/chromosome{chromosome}_haplotype{haplotype}_reference.fasta",
    shell:
        "graph_read_simulator prepare_simulation --chromosome {wildcards.chromosome} --haplotype {wildcards.haplotype} "
        "--vcf {input.vcf} --reference {input.reference} -o data/{wildcards.genome_build}/{wildcards.individual}/ "


def get_input_files_for_haplotype_sequence(wildcards):
    chromosomes = config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.size]["chromosomes"].split(",")
    return [
        f"data/{wildcards.genome_build}/{wildcards.individual}/chromosome{chromosome}_haplotype{wildcards.haplotype}_reference.fasta"
        for chromosome in chromosomes
    ]

rule make_haplotype_sequence_for_simulation:
    input:
        get_input_files_for_haplotype_sequence
    output:
        "data/{genome_build}/{individual}/{size}/haplotype{haplotype}.fa",
    shell:
        """
        cat {input} > {output}
        """

def get_haplotype_coordinate_maps(wildcards):
    chromosomes = config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.size][
        "chromosomes"].split(",")
    return [
        f"data/{wildcards.genome_build}/{wildcards.individual}/coordinate_map_chromosome{chromosome}_haplotype{wildcards.haplotype}.npz"
        for chromosome in chromosomes
    ]


rule merge_haplotype_coordinate_maps:
    input:
        get_haplotype_coordinate_maps
    output:
        "data/{genome_build}/{individual}/{size}/coordinate_map_haplotype{haplotype}.npz",
    run:
        from shared_memory_wrapper import to_file
        from graph_read_simulator.simulation import CoordinateMap, MultiChromosomeCoordinateMap
        chromosomes = config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.size]["chromosomes"].split(",")
        data = {}
        for chromosome, file in zip(chromosomes, input):
            data[chromosome] = CoordinateMap.from_file(file)
        to_file(MultiChromosomeCoordinateMap(data), output[0])


def get_mason_error_parameters(wildcards):
    profile = config["illumina_error_profiles"][wildcards.error_profile]
    return " --illumina-prob-deletion " + str(profile["deletion_prob"]) + \
            " --illumina-prob-insert " + str(profile["insertion_prob"]) + \
            " --illumina-prob-mismatch-scale " + str(profile["mismatch_scale"])

rule simulate_reads_for_chromosome_and_haplotype:
    input:
        haplotype_reference="{individual}/haplotype{haplotype}.fa"
    output:
        multiext("{individual}/whole_genome_single_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}", ".fq.gz", ".haplotype_truth.sam")
    conda:
        "../envs/mason.yml"
    params:
        error_parameters=get_mason_error_parameters,
        n_reads=lambda wildcards: int(wildcards.n_reads) // 2,
        mean_fragment_size=lambda wildcards: int(wildcards.read_length) * 3,
        min_fragment_size= lambda wildcards: int(wildcards.read_length) // 2,
        max_fragment_size= lambda wildcards: int(wildcards.read_length) * 6,
    threads:
        6
    shell:
        "mason_simulator -ir {input.haplotype_reference} -n {params.n_reads} -o {output[0]} -oa {output[1]} --num-threads 6 {params.error_parameters} "
        "--illumina-read-length {wildcards.read_length} "
        "--fragment-mean-size {params.mean_fragment_size} "
        "--fragment-min-size {params.min_fragment_size} "
        "--fragment-max-size {params.max_fragment_size} "


rule simulate_reads_for_chromosome_and_haplotype_paired_end:
    input:
        haplotype_reference="{individual}/haplotype{haplotype}.fa"
    output:
        reads1="{individual}/whole_genome_paired_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}-1.fq.gz",
        reads2="{individual}/whole_genome_paired_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}-2.fq.gz",
        truth1="{individual}/whole_genome_paired_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}.haplotype_truth.sam",
    conda:
        "../envs/mason.yml"
    params:
        error_parameters=get_mason_error_parameters,
        n_reads=lambda wildcards: int(wildcards.n_reads) // 4,  # divide by 4 for paird end since mason simulates n fragments
        mean_fragment_size= lambda wildcards: int(wildcards.read_length) * 3,
        min_fragment_size= lambda wildcards: int(wildcards.read_length) // 2,
        max_fragment_size= lambda wildcards: int(wildcards.read_length) * 6,
    threads:
        6
    shell:
        "mason_simulator -ir {input.haplotype_reference} -n {params.n_reads} -o {output.reads1} -or {output.reads2} -oa {output.truth1} --num-threads 6 {params.error_parameters} "
        "--illumina-read-length {wildcards.read_length} "
        "--fragment-mean-size {params.mean_fragment_size} "
        "--fragment-min-size {params.min_fragment_size} "
        "--fragment-max-size {params.max_fragment_size} "

rule merge_paired_end_reads:
    input:
        r1="{data}/{haplotype}-1.fq.gz",
        r2= "{data}/{haplotype}-2.fq.gz"
    output:
        merged="{data}/{haplotype,\d+}.fq.gz"
    params:
        compress_lvl=9,
    threads: 4
    wrapper:
        "v1.21.4/bio/seqtk/mergepe"


rule deinterleave_fastq:
    input:
        "{data}/reads.fq.gz"
    output:
        "{data}/reads{n}.fq.gz"
    params:
        extra="-{n}",
    conda:
        "../envs/seqtk.yml"
    shell:
        "seqtk seq -{wildcards.n} {input} | gzip -c > {output}"


# finds out whether each truth alignment covers a variant and adds that information
rule add_variant_info_to_truth_sam:
    input:
        truth_positions="{individual}/whole_genome_{pair}_end/{config}/{haplotype,\d+}.haplotype_truth.sam",
        coordinate_map="{individual}/coordinate_map_haplotype{haplotype}.npz",
    output:
        "{individual}/whole_genome_{pair}_end/{config}/{haplotype,\d+}.haplotype_truth.with_variant_info.sam",
    run:
        from shared_memory_wrapper import from_file
        coordinate_map = from_file(input.coordinate_map)
        out_file = open(output[0], "w")
        for line in open(input.truth_positions):
            if line.startswith("@"):
                out_file.write(line)
                continue

            l = line.split()
            chromosome = l[2]
            start = int(l[3])
            length = len(l[9])
            end = start + length
            n_variants = 0
            if coordinate_map.haplotype_has_variant_between(chromosome, start, end):
                n_variants = 1

            line = line.strip() + "\tNVARIANTS:i:" + str(n_variants) + "\n"
            out_file.write(line)


rule change_truth_alignments_to_reference_coordinates:
    input:
        truth_positions = "{individual}/whole_genome_{pair}_end/{config}/{haplotype,\d+}.haplotype_truth.with_variant_info.sam",
        coordinate_map="{individual}/coordinate_map_haplotype{haplotype}.npz",
    output:
        "{individual}/whole_genome_{pair}_end/{config}/{haplotype,\d+}.reference_coordinates.sam",
        #"data/simulated_reads/{dataset}/simulated_reads_haplotype{haplotype,\d+}.reference_coordinates.sam"
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
        haplotype0="{config}/0.fq.gz",
        haplotype1="{config}/1.fq.gz",
    output:
        "{config}/reads.fq.gz"
    shell:
        "zcat {input} | python scripts/assign_ids_to_fq.py | gzip > {output} "


rule merge_truth_alignments_for_haplotypes:
    input:
        haplotype0 = "{config}/0.reference_coordinates.sam",
        haplotype1 = "{config}/1.reference_coordinates.sam",
    output:
        "{config}/truth.sam"
    shell:
        "cat {input} | python scripts/assign_ids_to_sam.py  > {output} "


rule store_truth_alignments:
    input:
        sam="{data}/{n_reads}/truth.sam",
        #coordinate_map="data/simulated_reads/{dataset}/"
    output:
        "{data}/{n_reads,\d+}/truth.npz",
    shell:
        """
        cat {input.sam} | numpy_alignments store sam {output} {wildcards.n_reads}
        """

