configfile: "config/config.yaml"
configfile: "config/plots.yaml"

from mapping_benchmarking.config import *
from mapping_benchmarking.util import get_input_reads

wildcard_constraints:
    haplotype="0|1",


include: "rules/bowtie2.smk"
include: "rules/strobealign.smk"
include: "rules/minimap2.smk"
include: "rules/reference_genome.smk"
include: "rules/mason.smk"
include: "rules/read_simulation.smk"
include: "rules/chip_seq_simulation.smk"
include: "rules/chip_seq.smk"
include: "rules/analysis.smk"
include: "rules/reports.smk"
include: "rules/variant_calling.smk"
include: "rules/plotting.smk"
include: "rules/tests.smk"
include: "rules/bwa.smk"

