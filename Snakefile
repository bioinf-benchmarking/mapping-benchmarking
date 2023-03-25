configfile: "config/config.yaml"
configfile: "config/plots.yaml"
workflow.use_conda = True
from mapping_benchmarking.config import *


wildcard_constraints:
    n_threads="\d+",
    min_mapq="\d+",
    variant_filter="all|variants|nonvariants",
    genome_build="\w+",
    individual="\w+",
    dataset_size="small|medium|big",
    haplotype="0|1",
    peak_read_coverage="\d+"



#def paired_end_reads(wildcards=None):
#   return [f"data/{reference_genome}/{{config}}/reads" + n + ".fq.gz" for n in
#            ("1", "2")]


#def single_end_reads(wildcards=None):
#    return f"data/{reference_genome}/{{config}}/reads.fq.gz"


#def get_input_reads(wildcards):
#    return paired_end_reads() if "paired_end" in ''.join(wildcards) else single_end_reads()

def get_input_reads(wildcards):
    file = "reads.fq.gz" if not "paired_end" in wildcards.read_config else ["reads1.fq.gz", "reads2.fq.gz"]
    return GenericReads.path(file=file)


class Parameters:
    """
    Utility class for generating a snakemake string of all parameters
    Makes it possible to use {{parameters}} in a rule instead of manually
    specifying all parameters
    """
    def __init__(self, parameters, prefix=""):
        self.parameters = parameters
        self.prefix = prefix

    def __str__(self):
        splitter = "}/{"
        out = "{" + splitter.join(self.parameters) + "}"
        return self.prefix + out

    def __getitem__(self, item):
        return Parameters(self.parameters[item], prefix=self.prefix)

    def until(self, parameter):
        """ Gives all parameters until and including the given parameter"""
        assert parameter in self.parameters, "Parameter %s not in %s" % (parameter, self.parameters)
        return self[0:self.parameters.index(parameter)+1]

    def __call__(self, **filters):
        """
        Returns a format-string with parameters given in filters set to values instead of wildcards
        """
        out = []
        for parameter in self.parameters:
            if parameter in filters:
                out.append(str(filters[parameter]))
            else:
                out.append("{" + parameter + "}")
        return self.prefix + "/".join(out)

parameters_wgs = config["parameter_types_reference_genome"] + config["parameter_types_whole_genome_sequencing"]
parameters_chip_seq = config["parameter_types_reference_genome"] + config["parameter_types_chip_seq"]

parameters = Parameters(config["parameter_types"])

print("Parameters: ", type(parameters))
reference_genome = Parameters(config["parameter_types_reference_genome"])
wgs = Parameters(config["parameter_types_whole_genome_sequencing"])
chip_seq = Parameters(config["parameter_types_chip_seq"])


include: "rules/bowtie2.smk"
include: "rules/strobealign.smk"
include: "rules/minimap2.smk"
include: "rules/reference_genome.smk"
include: "rules/read_simulation.smk"
include: "rules/chip_seq_simulation.smk"
include: "rules/chip_seq.smk"
include: "rules/analysis.smk"
include: "rules/reports.smk"
include: "rules/variant_calling.smk"
include: "rules/plotting.smk"
include: "rules/tests.smk"
include: "rules/bwa.smk"



