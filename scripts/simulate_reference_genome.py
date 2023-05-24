import logging
import sys
import numpy as np
import typer
import bionumpy as bnp


def simulate_reference_genome(genome_size: int, n_chromosomes: int, out_file: str):
    sequences = bnp.simulate.simulate_sequences(
        "ACGT", {f"chr{i}": genome_size // n_chromosomes for i in range(n_chromosomes)})

    with bnp.open(out_file, "w") as f:
        f.write(sequences)


if __name__ == "__main__":
    typer.run(simulate_reference_genome)
