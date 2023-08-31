from mapping_benchmarking.config import ReferenceGenome

rule download_reference:
    output:
        "data/{reference}/reference.2bit"
    shell:
        "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.reference}/bigZips/{wildcards.reference}.2bit"


rule convert_reference_genome_to_fasta:
    input:
        "data/{genome_build}/reference.2bit"
    output:
        "data/{genome_build, ((?!simulated).)*}/reference.fa"
    wrapper:
        "v1.21.2/bio/ucsc/twoBitToFa"


rule simulate_reference_genome:
    output:
        "data/{genome_build, (simulated).+}/reference.fa"
    params:
        genome_size = lambda wildcards: config["genomes"][wildcards.genome_build]["genome_size"],
        n_chromosomes = lambda wildcards: config["genomes"][wildcards.genome_build]["n_chromosomes"],
    run:
        import bionumpy as bnp
        sequences = bnp.simulate.simulate_sequences(
            "ACGT",{f"chr{i}": params.genome_size // params.n_chromosomes for i in range(params.n_chromosomes)})

        with bnp.open(output[0],"w") as f:
            f.write(sequences)

    #"python3 scripts/simulate_reference_genome.py {params.genome_size} {params.n_chromosomes} {output}"



rule get_reference_genome_chromosome:
    input:
        ref="data/{genome_build}/reference.fa",
        index="data/{genome_build}/reference.fa.fai",
    output:
        "data/{genome_build}/reference_genome_by_chromosome/{chromosome}.fa.gz"
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools faidx {input} {wildcards.chromosome} | gzip -c > {output}"


rule samtools_index:
    input:
        "{sample}.fa",
    output:
        "{sample}.fa.fai",
    wrapper:
        "v2.6.0/bio/samtools/faidx"

rule samtools_index2:
    input:
        "{sample}.fasta",
    output:
        "{sample}.fasta.fai",
    wrapper:
        "v2.6.0/bio/samtools/faidx"


rule get_dataset_reference:
    input:
        ref = "data/{genome_build}/reference.fa",
        index = "data/{genome_build}/reference.fa.fai",
    output:
        ReferenceGenome.path()
    conda:
        "../envs/samtools.yml"
    params:
        regions=lambda wildcards: config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.dataset_size]["chromosomes"].replace(",", " ")
    shell:
        "samtools faidx {input.ref} {params.regions} > {output}"


rule tabix:
    input:
        "{prefix}.vcf.gz",
    output:
        "{prefix}.vcf.gz.tbi",
    params:
        "-p vcf",
    wrapper:
        "v1.21.6/bio/tabix/index"

