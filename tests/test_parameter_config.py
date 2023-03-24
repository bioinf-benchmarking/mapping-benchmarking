from mapping_benchmarking import config
from mapping_benchmarking.config import MappingF1Score, WholeGenomeReads, MappingAccuracy, GenericMappedReads, PeakCallingAccuracy
from snakemake.io import multiext


def test():
    print(PeakCallingAccuracy.from_flat_params().file_path())
    print(GenericMappedReads.as_output(method='bwa'))
    print(MappingAccuracy.path())
    MappingF1Score.from_flat_params(genome_build="hg38", genome="hg002", read_type="paired_end",
                                    read_length=150, n_reads=200, method="bwa").file_path()
    #print(WholeGenomeReads.from_flat_params().file_path())
    #print(MappingF1Score.from_flat_params(genome_build="hg38", genome="hg002", read_type="whole_genome_paired_end",
    #                                read_length=150, n_reads=200, method="bwa").file_path())
    #print(PairedEndReads.path())
    #print(PairedEndReads.path(file=[
    #    "reads1.fq.gz", "reads2.fq.gz"]))
    #print(AccuracyResult.path())
    #print(multiext(WholeGenomeSingleEnd.path(file_ending="") + "/{haplotype}", ".fq.gz", ".sam"))
    #print(parameter_config.WholeGenomeSingleEnd.as_output())
    #print(parameter_config.MappedReads.as_output())
    #print(parameter_config.ReferenceGenome.as_output())
    #print(ReferenceGenomeFile.path(file=["reference.fa", "reference.fa.fai"]))
