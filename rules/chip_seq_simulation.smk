import numpy as np
from mapping_benchmarking.config import ReferenceGenome, ChipSeq, SimulatedChipSeqPeaks

peak_size = 200


rule simulate_peaks:
    input:
        genome = ReferenceGenome.path(file_ending=".fa.fai")
    output:
        peaks = SimulatedChipSeqPeaks.path()
    script: "../scripts/simulate_peaks.py"


rule change_peak_coordinates_to_haplotype_coordinates:
    input:
        peaks = SimulatedChipSeqPeaks.path(),
        coordinate_map = ReferenceGenome.path(file_ending="") + "/coordinate_map_haplotype{haplotype}.npz"
    output:
        peaks = SimulatedChipSeqPeaks.path(file_ending="") + "_haplotype{haplotype}.bed",
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
        reference = ReferenceGenome.path(file_ending="") + "/haplotype{haplotype}.fa",
        reference_fai = ReferenceGenome.path(file_ending="") + "/haplotype{haplotype}.fa.fai",
        peaks = SimulatedChipSeqPeaks.path(file_ending="") + "_haplotype{haplotype}.bed",
    output:
        reads = ChipSeq.path() + "/{haplotype}.fq.gz"
    params:
        out_base_name = lambda wildcards, input, output: ".".join(output.reads.split(".")[:-1]),
        frac = peak_frac,
        n_reads=chip_seq_n_reads
    #conda:
    #    "../envs/chips.yml"
    shell:
        "chips simreads --spot {wildcards.spot} --numcopies 1000 --frac "
        "{params.frac} --seed 1 -t bed -c 5 -p {input.peaks} -f {input.reference} "
        "-o {params.out_base_name} --numreads {params.n_reads} --readlen {wildcards.read_length} "
        "&& gzip -c {params.out_base_name}.fastq > {output}"
