import numpy as np


rule simulate_peaks:
    input:
        genome=f"data/{reference_genome}/reference.fa.fai"
    output:
        peaks=f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks.bed"

    run:
        import bionumpy as bnp
        from bionumpy.simulate.intervals import simulate_fixed_size_uniform_intervals
        from bionumpy import Genome

        genome = Genome.from_file(input.genome)
        n_peaks = int(wildcards.n_peaks)
        simulated_peaks = simulate_fixed_size_uniform_intervals(genome, n_peaks, 500)

        # simulate scores
        scores = np.random.randint(0, 200, n_peaks)
        names = [f"peak{i}" for i in range(n_peaks)]
        strands = ["." for _ in range(n_peaks)]

        simulated_peaks = bnp.datatypes.Bed6(simulated_peaks.chromosome, simulated_peaks.start, simulated_peaks.stop,
            names, scores, strands)

        out_file = bnp.open(output.peaks, "w", buffer_type=bnp.Bed6Buffer)
        out_file.write(simulated_peaks)




rule change_peak_coordinates_to_haplotype_coordinates:
    input:
        peaks=f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks.bed",
        coordinate_map=f"data/{reference_genome}/coordinate_map_haplotype{{haplotype}}.npz"
    output:
        peaks = f"data/{reference_genome}/{chip_seq.until('n_peaks')}/simulated_peaks_haplotype{{haplotype,0|1}}.bed"
    shell:
        """
        graph_read_simulator liftover -i {input.peaks} -c {input.coordinate_map} -o {output.peaks}
        """


rule simulate_peak_reads:
    input:
        reference=f"data/{reference_genome}/reference.fa",
        peaks = f"data/{reference_genome}/{chip_seq.until('n_peaks')(read_type='chip_seq')}/simulated_peaks_haplotype{{haplotype}}.bed"
    output:
        reads=f"data/{reference_genome}/{chip_seq.until('n_peaks')(read_type='chip_seq')}/{{haplotype,0|1}}.fq"
    params:
        out_base_name = lambda wildcards, input, output: output.reads.split(".")[0],
        n_reads = lambda wildcards: int(wildcards.n_peaks) * 500
    conda:
        "../envs/chips.yml"
    shell:
        "chips simreads --seed 1 -t bed -c 5 -p {input.peaks} -f {input.reference} -o {params.out_base_name} --numreads {params.n_reads} --readlen {wildcards.read_length} && "
        "mv {params.out_base_name}.fastq {output}"

