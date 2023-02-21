
configfile: "config/config.yaml"
configfile: "config/plots.yaml"
workflow.use_conda = True


wildcard_constraints:
    n_threads="\d+",
    min_mapq="\d+",
    variant_filter="all|variants|nonvariants"


def paired_end_reads(wildcards=None):
    return [f"data/{parameters.until('n_reads')(read_type='whole_genome_paired_end')}/reads" + n + ".fq.gz" for n in
            ("1", "2")]


def single_end_reads(wildcards=None):
    return f"data/{parameters.until('n_reads')(read_type='whole_genome_single_end')}/reads.fq.gz"


def get_input_reads(wildcards):
    if "paired_end" in wildcards.read_type:
        return paired_end_reads()
    else:
        return single_end_reads()


class Parameters:
    """
    Utility class for generating a snakemake string of all parameters
    Makes it possible to use {{parameters}} in a rule instead of manually
    specifying all parameters
    """
    def __init__(self, parameters):
        self.parameters = parameters

    def __str__(self):
        splitter = "}/{"
        out = "{" + splitter.join(self.parameters) + "}"
        return out

    def __getitem__(self, item):
        return Parameters(self.parameters[item])

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
        return "/".join(out)


parameters = Parameters(config["parameter_types"])
print(str(parameters))

include: "rules/bwa.smk"
include: "rules/bowtie2.smk"
include: "rules/strobealign.smk"
include: "rules/minimap2.smk"
include: "rules/reference_genome.smk"
include: "rules/read_simulation.smk"
include: "rules/analysis.smk"
include: "rules/reports.smk"
include: "rules/variant_calling.smk"
include: "rules/plotting.smk"


rule test:
    output:
        touch(f"{parameters(genome='test123')}/test.txt")
    run:
        print(wildcards)



