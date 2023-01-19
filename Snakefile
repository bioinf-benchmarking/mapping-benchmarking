configfile: "config.yml"
workflow.use_conda = True

include: "rules/bwa.smk"
include: "rules/strobealign.smk"
include: "rules/minimap2.smk"
include: "rules/reference_genome.smk"
include: "rules/read_simulation.smk"
include: "rules/analysis.smk"

