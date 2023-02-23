

rule test:
    input:
        "data/hg38/hg002/small/whole_genome_single_end/medium_error/150/200/bwa/4/0/all/snps/f1_score.txt"
    output:
        touch("test.txt")
    run:
        score = float(open(input[0]).read().strip())
        assert score >= 0.90
        print("Tests passed")

