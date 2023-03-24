from snakehelp import parameters, result
from typing import Union, Literal
import snakehelp
snakehelp.set_data_folder("data/")


@parameters
class ReadMapper:
    method: Literal["bwa", "bowtie2", "minimap", "strobealign"] = "bwa"
    n_threads: int = 4


@parameters
class ReferenceGenome:
    file_ending = ".fa"
    genome_build: str = "hg38"
    individual: str = "hg002"
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
class FilteredWholeGenomeMappedReads:
    mapped_reads: WholeGenomeMappedReads
    min_mapq: int = 0
    variant_filter: Literal["all", "variants", "nonvariants"] = "all"


@parameters
class VariantCallingAccuracy:
    mapped_reads: FilteredWholeGenomeMappedReads
    variant_type: Literal["snps", "indels", "all"]
    accuracy_type: Literal["recall", "one_minus_precision", "f1_score"]
    file_ending = ".txt"


@parameters
class MappingAccuracy:
    filtered_mapped_reads: FilteredWholeGenomeMappedReads
    accuracy_type: Literal["MappingRecall", "MappingOneMinusPrecision", "MappingF1Score"]
    file_ending = ".txt"


@result
class MappingRecall:
    filtered_mapped_reads: FilteredWholeGenomeMappedReads


@result
class MappingOneMinusPrecision:
    filtered_mapped_reads: FilteredWholeGenomeMappedReads


@result
class MappingF1Score:
    filtered_mapped_reads: FilteredWholeGenomeMappedReads


@result
class Runtime:
    filtered_mapped_reads: FilteredWholeGenomeMappedReads


@result
class MemoryUsage:
    filtered_mapped_reads: FilteredWholeGenomeMappedReads


@result
class PeakCallingAccuracy:
    called_peaks: CalledPeaks

