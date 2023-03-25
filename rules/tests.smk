from mapping_benchmarking.config import *


rule test_that_mappers_work:
    input:
        WholeGenomeMappedReads.from_flat_params(method="bowtie2", n_reads=1000).file_path(),
        WholeGenomeMappedReads.from_flat_params(method="minimap", n_reads=1000).file_path(),
        WholeGenomeMappedReads.from_flat_params(method="strobealign", n_reads=1000).file_path(),
    output:
        touch("test_mappers.txt")


rule test_runtime:
    input:
        Runtime.from_flat_params(method="bwa", n_reads=1000).file_path()
    output:
        touch("test_runtime.txt")


rule test_filtering_on_mapq:
    input:
        FilteredWholeGenomeMappedReads.from_flat_params(method="bwa", n_reads=1000, mapq=10).file_path()
    output:
        touch("test_filtering_on_mapq.txt")


rule test_variant_calling:
    input:
        VariantCalls.from_flat_params(method="bwa", mapq=10).file_path()
    output:
        touch("test_variant_calling.txt")


rule test_variant_calling_accuracy:
    input:
        VariantCallingRecall.from_flat_params().file_path()
    output:
        touch("test_variant_calling_accuracy.txt")



rule test_accuracy:
    input:
        #"data/hg38/hg002/small/whole_genome_single_end/medium_error/150/200/bwa/4/0/all/f1_score.txt"
        MappingF1Score.from_flat_params(genome_build="hg38",genome="hg002",read_type="single_end",read_length=150,n_reads=200,method="bwa").file_path()
    output:
        touch("test_accuracy.txt")
    run:
        score = float(open(input[0]).read().strip())
        assert score >= 0.90
        print("Tests passed")


rule test_accuracy_paired_end:
    input:
        MappingF1Score.from_flat_params(genome_build="hg38", genome="hg002", read_type="paired_end", read_length=150, n_reads=200, method="bwa").file_path()
    output:
        touch("test_accuracy_paired_end.txt")
    run:
        score = float(open(input[0]).read().strip())
        assert score >= 0.95
        print("Tests passed")



rule test_peak_calling_accuracy:
    input:
        PeakCallingAccuracy.from_flat_params(genome_build="sacCer3", individual="simulated").file_path()
        #"data/sacCer3/simulated/small/chip_seq/medium_error/75/0.3/40/50/bwa/4/accuracy.txt"
    output:
        touch("test_peak_calling_accuracy.txt")
    run:
        score = float(open(input[0]).read().strip())
        assert score >= 0.90


rule test:
    input:
        "test_accuracy.txt",
        "test_peak_calling_accuracy.txt",
    output:
        touch("test.txt")