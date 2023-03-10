

rule test_accuracy:
    input:
        "data/hg38/hg002/small/whole_genome_single_end/medium_error/150/200/bwa/4/0/all/snps/f1_score.txt"
    output:
        touch("test.txt")
    run:
        score = float(open(input[0]).read().strip())
        assert score >= 0.90
        print("Tests passed")



rule test:
    input:
        "data/sacCer3/simulated/small/whole_genome_single_end/medium_error/150/10000/bwa/4/mapped.npz",
        "data/sacCer3/simulated/small/whole_genome_single_end/medium_error/150/10000/bwa/4/0/all/snps/f1_score.txt",
        "data/sacCer3/simulated/small/chip_seq/medium_error/1000/1.fq.gz"
    output:
        touch("test2.txt")