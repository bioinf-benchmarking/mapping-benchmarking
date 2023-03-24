from mapping_benchmarking.config import MappingF1Score


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