
configfile: "config/config.yaml"
configfile: "config/plots.yaml"
workflow.use_conda = True

class Parameters:
    """
    Utility class for generating a snakemake string of all parameters
    Makes it possible to use {{parameters}} in a rule instead of manually
    specifying all parameters
    """
    def __init__(self, parameters):
        self.parameters = parameters

    def __str__(self):
        return "{" + "}/{".join(self.parameters) + "}"

    def __getitem__(self, item):
        return Parameters(self.parameters[item])

    def until(self, parameter):
        """ Gives all parameters until and including the given parameter"""
        assert parameter in self.parameters
        return self[0:self.parameters.index(parameter)+1]

    def __call__(self, **filters):
        """
        Returns a format-string with parameters given in filters set to values instead of wildcards
        """
        out = []
        for parameter in self.parameters:
            if parameter in filters:
                out.append(filters[parameter])
            else:
                out.append("{" + parameter + "}")
        return "/".join(out)


parameters = Parameters(config["parameter_types"])

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



