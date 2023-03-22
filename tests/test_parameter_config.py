from mapping_benchmarking import parameter_config
from mapping_benchmarking.parameter_config import WholeGenomeSingleEnd
from snakemake.io import multiext


def test():
    print(multiext(WholeGenomeSingleEnd.path(file_ending="") + "/{haplotype}", ".fq.gz", ".sam"))
    #print(parameter_config.WholeGenomeSingleEnd.as_output())
    #print(parameter_config.MappedReads.as_output())
    #print(parameter_config.ReferenceGenome.as_output())
    #print(ReferenceGenomeFile.path(file=["reference.fa", "reference.fa.fai"]))
