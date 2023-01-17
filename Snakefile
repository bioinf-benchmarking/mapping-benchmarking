configfile: "config.yml"


include: "rules/bwa.smk"
include: "rules/reference_genome.smk"
include: "rules/read_simulation.smk"
include: "rules/analysis.smk"
