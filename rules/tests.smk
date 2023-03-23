

rule test_accuracy:
    input:
        "data/hg38/hg002/small/whole_genome_single_end/medium_error/150/200/bwa/4/0/all/snps/f1_score.txt"
    output:
        touch("test_accuracy.txt")
    run:
        score = float(open(input[0]).read().strip())
        assert score >= 0.90
        print("Tests passed")



rule test_peak_calling_accuracy:
    input:
        "data/sacCer3/simulated/small/chip_seq/medium_error/75/0.3/40/50/bwa/4/accuracy.txt"
    output:
        touch("test_peak_calling_accuracy.txt")
    run:
        score = float(open(input[0]).read().strip())
        assert score >= 0.90


rule test:
    input:
        #"test_accuracy.txt",
        "test_peak_calling_accuracy.txt",
        "data/sacCer3/simulated/small/whole_genome_single_end/medium_error/150/10000/bwa/4.npz",
        "data/sacCer3/simulated/small/whole_genome_single_end/medium_error/150/10000/bwa/4/0/all/snps/f1_score.txt",
    output:
        touch("test.txt")