

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
        haplotype_reference="data/{genome_build}/{individual}/chromosome{chromosome}_haplotype{haplotype}_reference.fasta"
    shell:
        "graph_read_simulator prepare_simulation --chromosome {wildcards.chromosome} --haplotype {wildcards.haplotype} "
        "--vcf {input.vcf} --reference {input.reference} -o data/{wildcards.genome_build}/{wildcards.individual}/"


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
            print("Chromosome %s, %s" % (chromosome, file))
            data[chromosome] = CoordinateMap.from_file(file)
        to_file(MultiChromosomeCoordinateMap(data), output[0])


def get_mason_error_parameters(wildcards):
    profile = config["illumina_error_profiles"][wildcards.error_profile]
    return " --illumina-prob-deletion " + str(profile["deletion_prob"]) + \
            " --illumina-prob-insert " + str(profile["insertion_prob"]) + \
            " --illumina-prob-mismatch-scale " + str(profile["mismatch_scale"])

def get_art_error_parameters(wildcards):
    profile = config["illumina_error_profiles"][wildcards.error_profile]
    return \
            " --delRate " + str(profile["deletion_prob"]) + \
            " --delRate2 " + str(profile["deletion_prob"]) + \
            " --insRate " + str(profile["insertion_prob"]) + \
            " --insRate2 " + str(profile["insertion_prob"]) + \
            " -qs " + str(1 / profile["mismatch_scale"]) + \
            " -qs2 " + str(1 / profile["mismatch_scale"])


"""
rule simulate_reads_for_chromosome_and_haplotype:
    input:
        haplotype_reference="{individual}/haplotype{haplotype}.fa"
    output:
        multiext("{individual}/whole_genome_single_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}", ".fq", ".sam")
    conda:
        "../envs/mason.yml"
    params:
        error_parameters=get_mason_error_parameters,
        n_reads=lambda wildcards: int(wildcards.n_reads) // 2,
        mean_fragment_size=lambda wildcards: int(wildcards.read_length) * 3,
        min_fragment_size= lambda wildcards: int(wildcards.read_length) // 2,
        max_fragment_size= lambda wildcards: int(wildcards.read_length) * 6,
    threads:
        2
    shell:
        "mason_simulator -ir {input.haplotype_reference} -n {params.n_reads} -o {output[0]} -oa {output[1]} --num-threads 2 {params.error_parameters} "
        "--illumina-read-length {wildcards.read_length} "
        "--fragment-mean-size {params.mean_fragment_size} "
        "--fragment-min-size {params.min_fragment_size} "
        "--fragment-max-size {params.max_fragment_size} "
"""

def get_genome_size(wildcards, input, output):
    with open(input.haplotype_reference_fai) as f:
        size = sum((int(l.strip().split()[1]) for l in f))
    print("Genome size:", size)
    return size


def get_coverage(wildcards, input, output):
    genome_size = config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.dataset_size]["genome_size"]  # get_genome_size(wildcards, input, output)
    coverage = int(wildcards.n_reads) * int(wildcards.read_length) / genome_size
    return coverage / 2


rule simulate_reads_for_chromosome_and_haplotype_art:
    input:
        haplotype_reference=f"data/{parameters.until('dataset_size')}/haplotype{{haplotype}}.fa",
        #fai=f"data/{parameters.until('dataset_size')}/haplotype{{haplotype}}.fa.fai",
        #haplotype_reference_fai="{individual}/haplotype{haplotype}.fa.fai",
    output:
        multiext(f"data/{parameters.until('n_reads')(read_type='whole_genome_single_end')}/{{haplotype,\d+}}", ".fq", ".sam")
    conda:
        "../envs/art.yml"
    params:
        error_parameters=get_art_error_parameters,
        coverage=get_coverage,
        mean_fragment_size=lambda wildcards: int(wildcards.read_length) * 3,
        min_fragment_size= lambda wildcards: int(wildcards.read_length) // 2,
        max_fragment_size= lambda wildcards: int(wildcards.read_length) * 6,
        output_base_name=lambda wildcards, input, output: output[0].split(".")[0]
    threads:
        2
    shell:
        "art_illumina -ss MSv3 -sam -i {input.haplotype_reference} -f {params.coverage} -o {params.output_base_name} -l {wildcards.read_length} {params.error_parameters} "


rule simulate_reads_for_chromosome_and_haplotype_paired_end_art:
    input:
        haplotype_reference=f"data/{parameters.until('dataset_size')}/haplotype{{haplotype}}.fa",
        #fai=f"data/{parameters.until('dataset_size')}/haplotype{{haplotype}}.fa.fai",
        #haplotype_reference="{individual}/haplotype{haplotype}.fa",
        #haplotype_reference_fai="{individual}/haplotype{haplotype}.fa.fai",
    output:
        multiext(f"data/{parameters.until('n_reads')(read_type='whole_genome_paired_end')}/{{haplotype,\d+}}", "-1.fq", "-2.fq", "-.sam")
        #multiext("{individual}/whole_genome_paired_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}", "-1.fq", "-2.fq", "-.sam")
    conda:
        "../envs/art.yml"
    params:
        error_parameters=get_art_error_parameters,
        coverage=get_coverage,
        mean_fragment_size=lambda wildcards: int(wildcards.read_length) * 3,
        min_fragment_size= lambda wildcards: int(wildcards.read_length) // 2,
        max_fragment_size= lambda wildcards: int(wildcards.read_length) * 6,
        std_fragment_size= lambda wildcards: int(wildcards.read_length) // 10 ,
        output_base_name=lambda wildcards, input, output: output[0].split(".")[0][0:-1]
    threads:
        2
    shell:
        "art_illumina -ss MSv3 -p -sam -i {input.haplotype_reference} -f {params.coverage} -o {params.output_base_name} "
        "-l {wildcards.read_length} -m {params.mean_fragment_size} -s {params.std_fragment_size} {params.error_parameters} "


# hack to get paired end rule to give same as single end
rule fix_sam_file_name:
    input:
        f"data/{parameters.until('n_reads')(read_type='whole_genome_paired_end')}/{{haplotype,\d+}}-.sam"
        #"{dir}/{haplotype,\d+}-.sam"
    output:
        f"data/{parameters.until('n_reads')(read_type='whole_genome_paired_end')}/{{haplotype,\d+}}.sam"
        #"{dir}/{haplotype,\d+}.sam"
    shell: "cp {input} {output}"


rule gz:
    input: "{file}.fq"
    output: "{file}.fq.gz"
    shell: "gzip -c {input} > {output}"


"""
rule simulate_reads_for_chromosome_and_haplotype_paired_end:
    input:
        haplotype_reference="{individual}/haplotype{haplotype}.fa"
    output:
        reads1=temp("{individual}/whole_genome_paired_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}-1.fq.gz"),
        reads2=temp("{individual}/whole_genome_paired_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}-2.fq.gz"),
        truth1=temp("{individual}/whole_genome_paired_end/{error_profile}/{read_length}/{n_reads}/{haplotype,\d+}.haplotype_truth.sam"),
    conda:
        "../envs/mason.yml"
    params:
        error_parameters=get_mason_error_parameters,
        n_reads=lambda wildcards: int(wildcards.n_reads) // 4,  # divide by 4 for paird end since mason simulates n fragments
        mean_fragment_size= lambda wildcards: int(wildcards.read_length) * 3,
        min_fragment_size= lambda wildcards: int(wildcards.read_length) // 2,
        max_fragment_size= lambda wildcards: int(wildcards.read_length) * 6,
    threads:
        2
    shell:
        "mason_simulator -ir {input.haplotype_reference} -n {params.n_reads} -o {output.reads1} -or {output.reads2} -oa {output.truth1} --num-threads 2 {params.error_parameters} "
        "--illumina-read-length {wildcards.read_length} "
        "--fragment-mean-size {params.mean_fragment_size} "
        "--fragment-min-size {params.min_fragment_size} "
        "--fragment-max-size {params.max_fragment_size} "
"""

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
        truth_positions="{individual}/whole_genome_{pair}_end/{config}/{haplotype,\d+}.sam",
        coordinate_map="{individual}/coordinate_map_haplotype{haplotype}.npz",
    output:
        "{individual}/whole_genome_{pair}_end/{config}/{haplotype,\d+}.haplotype_truth.with_variant_info.sam",
    script: "../scripts/add_variant_info_to_truth_sam.py"



rule change_truth_alignments_to_reference_coordinates:
    input:
        truth_positions = "{individual}/whole_genome_{pair}_end/{config}/{haplotype,\d+}.haplotype_truth.with_variant_info.sam",
        coordinate_map="{individual}/coordinate_map_haplotype{haplotype}.npz",
    output:
        "{individual}/whole_genome_{pair}_end/{config}/{haplotype,\d+}.reference_coordinates.sam"
        #"data/simulated_reads/{dataset}/simulated_reads_haplotype{haplotype,\d+}.reference_coordinates.sam"
    script: "../scripts/change_truth_alignments_to_reference_coordinates.py"



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


rule count_reads:
    input:
        f"data/{parameters.until('n_reads')}/truth.sam"
    output:
        f"data/{parameters.until('n_reads')}/n_reads.txt"
    shell:
        "wc -l {input} | cut -f 1 -d ' ' > {output}"


rule store_truth_alignments:
    input:
        reads="{data}/{n_reads}/reads.fq.gz",  # not necessary
        sam="{data}/{n_reads}/truth.sam",
        n_reads="{data}/{n_reads}/n_reads.txt",
    params:
        n_reads=lambda wildcards, input, output: open(input.n_reads).read().strip()
    output:
        "{data}/{n_reads,\d+}/truth.npz",
    shell:
        """
        cat {input.sam} | numpy_alignments store sam {output} {params.n_reads}
        """

