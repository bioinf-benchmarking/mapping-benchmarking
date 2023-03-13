import numpy as np
import bionumpy as bnp
from bionumpy.simulate.intervals import simulate_fixed_size_uniform_intervals
from bionumpy import Genome

peak_size = 100

genome = Genome.from_file(snakemake.input.genome)
n_peaks = int(snakemake.wildcards.n_peaks)
print("Peak size: %d" % peak_size)
simulated_peaks = simulate_fixed_size_uniform_intervals(genome, n_peaks, peak_size)

# simulate scores
scores = np.random.randint(150, 200, n_peaks)
names = [f"peak{i}" for i in range(n_peaks)]
strands = ["." for _ in range(n_peaks)]

simulated_peaks = bnp.datatypes.Bed6(simulated_peaks.chromosome, simulated_peaks.start, simulated_peaks.stop,
    names, scores, strands)

out_file = bnp.open(snakemake.output.peaks, "w", buffer_type=bnp.Bed6Buffer)
out_file.write(simulated_peaks)
