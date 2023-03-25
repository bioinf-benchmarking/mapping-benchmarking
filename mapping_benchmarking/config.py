from snakehelp import parameters, result
from typing import Union, Literal
import snakehelp
snakehelp.set_data_folder("data/")


@parameters
class ReadMapper:
    method: Literal["bwa", "bowtie2", "minimap", "strobealign"] = "bwa"
    n_threads: int = 4


@parameters
class GenomeBuild:
    genome_build: str = "hg38"


@parameters
class Individual:
    genome_build: GenomeBuild
    individual: str = "hg002"


@parameters
class ReferenceGenome:
    file_ending = ".fa"
    individual: Individual
    dataset_size: Literal["small", "medium", "large"] = "small"


@parameters
class WholeGenomeConfig:
    read_type: Literal["single_end", "paired_end"] = "single_end"
    error_profile: str = "medium_error"
    read_length: int = 150
    n_reads: int = 1000


@parameters
class ChipSeqConfig:
    experiment_type: Literal["chip_seq"] = "chip_seq"
    error_profile: Literal["small_error", "medium_error", "high_error"] = "medium_error"
    read_length: int = 150
    spot: float = 0.2  # chips term: Ratio of reads that are in peaks
    peak_read_coverage: int = 50
    n_peaks: int = 200


@parameters
class ChipSeq:
    reference_genome: ReferenceGenome
    config: ChipSeqConfig


@parameters
class CalledPeaks:
    experiment: ChipSeq
    mapper: ReadMapper


@parameters
class SimulatedChipSeqPeaks:
    config: ChipSeq
    file: Literal["simulated_peaks"] = "simulated_peaks"
    file_ending = ".bed"



@parameters
class WholeGenomeReads:
    reference_genome: ReferenceGenome
    read_config: WholeGenomeConfig


@parameters
class PairedEndReads:
    reference_genome: ReferenceGenome
    read_config: WholeGenomeConfig
    file: Literal["reads1.fq.gz", "reads2.fq.gz"]


@parameters
class GenericReads:
    reference_genome: ReferenceGenome
    read_config: Union[WholeGenomeConfig, ChipSeqConfig]
    file: Literal["reads.fq.gz", "reads1.fq.gz", "reads2.fq.gz"]


@parameters
class WholeGenomeMappedReads:
    reads: WholeGenomeReads
    mapper: ReadMapper
    file_ending = ".bam"


@parameters
class GenericMappedReads:
    reference_genome: ReferenceGenome
    read_config: Union[WholeGenomeConfig, ChipSeqConfig]
    mapper: ReadMapper
    file_ending = ".bam"



@parameters
class MapQFilteredWholeGenomeMappedReads:
    mapped_reads: WholeGenomeMappedReads
    min_mapq: int = 0
    file_ending = ".bam"


@parameters
class FilteredWholeGenomeMappedReads:
    mapq_filtered_mapped_reads: MapQFilteredWholeGenomeMappedReads
    variant_filter: Literal["all", "variants", "nonvariants"] = "all"


@parameters
class VariantCalls:
    mapped_reads: MapQFilteredWholeGenomeMappedReads
    file_ending = ".vcf.gz"


@parameters
class FilteredVariantCalls:
    variant_calls: VariantCalls
    variant_calling_type: Literal["snps", "indels", "all"] = "all"


@parameters
class VariantCallingAccuracy:
    filtered_variant_calls: FilteredVariantCalls
    accuracy_type: Literal["VariantCallingRecall", "VariantCallingOneMinusPrecision", "VariantCallingF1Score"]
    file_ending = ".txt"


@result
class VariantCallingRecall:
    variant_calls: VariantCalls


@result
class VariantCallingF1Score:
    variant_calls: VariantCalls


@result
class VariantCallingOneMinusPrecision:
    variant_calls: VariantCalls


@parameters
class MappingAccuracy:
    filtered_mapped_reads: MapQFilteredWholeGenomeMappedReads
    accuracy_type: Literal["MappingRecall", "MappingOneMinusPrecision", "MappingF1Score"]
    file_ending = ".txt"


@result
class MappingRecall:
    filtered_mapped_reads: MapQFilteredWholeGenomeMappedReads


@result
class MappingOneMinusPrecision:
    filtered_mapped_reads: MapQFilteredWholeGenomeMappedReads


@result
class MappingF1Score:
    filtered_mapped_reads: MapQFilteredWholeGenomeMappedReads


@result
class Runtime:
    filtered_mapped_reads: MapQFilteredWholeGenomeMappedReads


@result
class MemoryUsage:
    filtered_mapped_reads: MapQFilteredWholeGenomeMappedReads


@result
class PeakCallingAccuracy:
    called_peaks: CalledPeaks

