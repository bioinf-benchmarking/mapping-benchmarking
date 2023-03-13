import numpy as np

peak_size = 200


rule simulate_peaks:
    input:
        genome=f"data/{reference_genome}/reference.fa.fai"
    output:
        peaks=f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks.bed"
    script: "../scripts/simulate_peaks.py"


rule make_chipulate_input_file:
    input:
        peaks = f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks.bed"
    output:
        tsv = f"data/{reference_genome}/{chip_seq.until('n_peaks')}/chipulate_input.tsv"
    run:
        import bionumpy as bnp
        peaks = bnp.open(input.peaks).read()
        with open(output.tsv, "w") as f:
            f.write("chr\tstart\tend\tp_ext\tp_amp\tenergy_A\n")
            for peak in peaks:
                p_ext = np.random.randint(40, 60) / 100
                p_amp = np.random.randint(10, 60) / 100
                energy_A = np.random.randint(10, 60) / 100
                f.write(f"{peak.chromosome}\t{peak.start}\t{peak.stop}\t{p_ext}\t{p_amp}\t{energy_A}\n")


rule change_peak_coordinates_to_haplotype_coordinates:
    input:
        peaks=f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks.bed",
        coordinate_map=f"data/{reference_genome}/coordinate_map_haplotype{{haplotype}}.npz"
    output:
        peaks = f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks_haplotype{{haplotype,0|1}}.bed"
    shell:
        """
        graph_read_simulator liftover -i {input.peaks} -c {input.coordinate_map} -o {output.peaks} --reverse True
        """


def genome_size(wildcards):
    return config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.dataset_size]["genome_size"]


def peak_regions_length(wildcards):
    return int(wildcards.n_peaks) * peak_size


def peak_frac(wildcards):
    genome_size = int(config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.dataset_size]["genome_size"])
    n_peaks = int(wildcards.n_peaks)
    peak_size = 200
    frac = n_peaks * peak_size / genome_size
    print("Fraction of genome covered by peaks: %.5f" % frac)
    return frac

def chip_seq_n_reads(wildcards):
    spot = float(wildcards.spot)
    coverage_in_peaks = int(wildcards.peak_read_coverage)
    # x0.5 since this is one out of two haplotypes
    n_reads_in_peaks = 0.5 * coverage_in_peaks * peak_regions_length(wildcards) / int(wildcards.read_length)
    n_reads_outside_peaks = n_reads_in_peaks / spot - n_reads_in_peaks

    print("N reads in peaks: %d" % n_reads_in_peaks)
    print("N reads outside peaks: %d" % n_reads_outside_peaks)
    total_reads = n_reads_in_peaks + n_reads_outside_peaks
    print("Total reads to simulate: %d" % total_reads)
    return int(total_reads)


rule simulate_peak_reads:
    input:
        reference=f"data/{reference_genome}/haplotype{{haplotype}}.fa",
        reference_fai=f"data/{reference_genome}/haplotype{{haplotype}}.fa.fai",
        peaks = f"data/{reference_genome}/{chip_seq.until('n_peaks')(read_type='chip_seq')}/simulated_peaks_haplotype{{haplotype}}.bed"
    output:
        reads=f"data/{reference_genome}/{chip_seq.until('n_peaks')(read_type='chip_seq')}/{{haplotype,0|1}}.fq.gz"
    params:
        out_base_name = lambda wildcards, input, output: ".".join(output.reads.split(".")[:-1]),
        frac = peak_frac,
        n_reads=chip_seq_n_reads
    conda:
        "../envs/chips.yml"
    shell:
        "chips simreads --spot {wildcards.spot} --numcopies 1000 --frac "
        "{params.frac} --seed 1 -t bed -c 5 -p {input.peaks} -f {input.reference} "
        "-o {params.out_base_name} --numreads {params.n_reads} --readlen {wildcards.read_length} "
        "&& gzip -c {params.out_base_name}.fastq > {output}"

"""
def isChip_ref_input(wildcards):
    chromosomes = config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.dataset_size]["chromosomes"].split(",")
    return ["data/" + wildcards.genome_build + "/reference_genome_by_chromosome/" + chromosome + ".fa.gz" for chromosome in chromosomes]


rule simulate_peak_reads_with_isChip:
    input:
        peaks = f"data/{reference_genome}/{chip_seq.until('n_peaks')(read_type='chip_seq')}/simulated_peaks_haplotype{{haplotype}}.bed",
        ref=isChip_ref_input
    output:
        reads=f"data/{reference_genome}/{chip_seq.until('n_peaks')(read_type='chip_seq')}/{{haplotype}}.fq"
    params:
        out_folder=lambda wildcards, input, output: "/".join(output.reads.split("/")[:-1]),
        ref_folder=lambda wildcards, input, output: "/".join(input[1].split("/")[:-1]),
    shell:
        "isChIP -O {params.out_folder} -g {params.ref_folder}/ -f fq,sam {input.peaks}"
"""