from mapping_benchmarking.config import ReferenceGenome, WholeGenomeReads, Individual, GenomeBuild



rule download_truth_regions:
    output:
        Individual.path() + "/truth_regions.bed"
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
        "graph_read_simulator prepare_simulation --chromosome '{wildcards.chromosome}' --haplotype {wildcards.haplotype} "
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
        haplotype_reference=ReferenceGenome.path(file_ending="") + "/haplotype{haplotype}.fa",
    output:
        multiext(WholeGenomeReads.path(read_type="single_end", file_ending="") + "/{haplotype}", ".fq.gz", ".sam")
    conda:
        "../envs/art.yml"
    params:
        error_parameters=get_art_error_parameters,
        coverage=get_coverage,
        mean_fragment_size=lambda wildcards: int(wildcards.read_length) * 3,
        min_fragment_size= lambda wildcards: int(wildcards.read_length) // 2,
        max_fragment_size= lambda wildcards: int(wildcards.read_length) * 6,
        output_base_name=lambda wildcards, input, output: output[0].split(".")[0],
        reads_base= lambda wildcards,input,output: output[0].replace(".gz","")
    threads:
        2
    shell:
        "art_illumina -ss MSv3 -sam -i {input.haplotype_reference} -f {params.coverage} -o {params.output_base_name} -l {wildcards.read_length} {params.error_parameters} "
        "&& gzip {params.reads_base}"


rule simulate_reads_for_chromosome_and_haplotype_paired_end_art:
    input:
        haplotype_reference=ReferenceGenome.path(file_ending="") + "/haplotype{haplotype}.fa",
    output:
        multiext(WholeGenomeReads.path(read_type="paired_end", file_ending="") + "/{haplotype}", "-1.fq.gz", "-2.fq.gz", "-.sam")
    conda:
        "../envs/art.yml"
    params:
        error_parameters=get_art_error_parameters,
        coverage=get_coverage,
        mean_fragment_size=lambda wildcards: int(wildcards.read_length) * 3,
        min_fragment_size= lambda wildcards: int(wildcards.read_length) // 2,
        max_fragment_size= lambda wildcards: int(wildcards.read_length) * 6,
        std_fragment_size= lambda wildcards: int(wildcards.read_length) // 10 ,
        output_base_name=lambda wildcards, input, output: output[0].split(".")[0][0:-1],
        reads1_base=lambda wildcards, input, output: output[0].replace(".gz", ""),
        reads2_base= lambda wildcards,input,output: output[1].replace(".gz", "")
    threads:
        2
    shell:
        "art_illumina -ss MSv3 -p -sam -i {input.haplotype_reference} -f {params.coverage} -o {params.output_base_name} "
        "-l {wildcards.read_length} -m {params.mean_fragment_size} -s {params.std_fragment_size} {params.error_parameters} && "
        "gzip {params.reads1_base} && gzip {params.reads2_base}"


# hack to get paired end rule to give same as single end
rule fix_sam_file_name:
    input:
        WholeGenomeReads.path(read_type="paired_end") + "/{haplotype}-.sam"
    output:
        WholeGenomeReads.path(read_type="paired_end") + "/{haplotype}.sam"

    shell: "cp {input} {output}"


#rule gz:
#    input: "{file}.fq"
#    output: "{file}.fq.gz"
#    shell: "gzip -c {input} > {output}"



rule merge_paired_end_reads:
    input:
        r1 = WholeGenomeReads.path(read_type="paired_end") + "/{haplotype}-1.fq.gz",
        r2 = WholeGenomeReads.path(read_type="paired_end") + "/{haplotype}-2.fq.gz",
    output:
        merged = WholeGenomeReads.path(read_type="paired_end") + "/{haplotype}.fq.gz",
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
        truth_positions=WholeGenomeReads.path(file_ending="") + "/{haplotype}.sam",
        coordinate_map=ReferenceGenome.path(file_ending="") + "/coordinate_map_haplotype{haplotype}.npz",
    output:
        sam = WholeGenomeReads.path(file_ending="") + "/{haplotype}.haplotype_truth.with_variant_info.sam",
        txt = WholeGenomeReads.path(file_ending="") + "/{haplotype}.n_variants.txt",
    script: "../scripts/add_variant_info_to_truth_sam.py"



rule change_truth_alignments_to_reference_coordinates:
    input:
        truth_positions=WholeGenomeReads.path(file_ending="") + "/{haplotype}.haplotype_truth.with_variant_info.sam",
        coordinate_map=ReferenceGenome.path(file_ending="") + "/coordinate_map_haplotype{haplotype}.npz",
    output:
        truth_positions=WholeGenomeReads.path(file_ending="") + "/{haplotype}.reference_coordinates.sam",
    script: "../scripts/change_truth_alignments_to_reference_coordinates.py"



rule merge_simulated_reads_for_haplotypes:
    input:
        haplotype0="{config}/0.fq.gz",
        haplotype1="{config}/1.fq.gz",
    output:
        "{config}/reads.fq.gz"
    shell:
        "zcat {input} | python scripts/assign_ids_to_fq.py | gzip -c > {output} "


rule merge_truth_alignments_for_haplotypes:
    input:
        haplotype0 = "{config}/0.reference_coordinates.sam",
        haplotype1 = "{config}/1.reference_coordinates.sam",
        haplotype0_n_variants = "{config}/0.n_variants.txt",
        haplotype1_n_variants= "{config}/1.n_variants.txt",
    output:
        bam="{config}/truth.bam",
        n_variants="{config}/n_variants.txt"
    conda:
        "../envs/samtools.yml"
    shell:
        "cat {input.haplotype0} {input.haplotype1} | python scripts/assign_ids_to_sam.py  | samtools view -o {output.bam} - && "
        "cat {input.haplotype0_n_variants} {input.haplotype1_n_variants} > {output.n_variants}"


rule store_truth_alignments:
    input:
        bam="{data}/truth.bam",
        n_variants="{data}/n_variants.txt",
    output:
        "{data}/truth.npz",
    shell:
        "numpy_alignments store bam -i {input.bam} -n {input.n_variants} {output} -1"

