from mapping_benchmarking.parameter_config import ReferenceGenome

rule download_reference:
    output:
        "data/{reference}/reference.2bit"
    shell:
        "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.reference}/bigZips/{wildcards.reference}.2bit"


rule convert_reference_genome_to_fasta:
    input:
        "data/{genome_build}/reference.2bit"
    output:
        "data/{genome_build}/reference.fa"
    wrapper:
        "v1.21.2/bio/ucsc/twoBitToFa"


rule get_reference_genome_chromosome:
    input:
        "data/{genome_build}/reference.fa"
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
        "v1.21.2/bio/samtools/faidx"

rule samtools_index2:
    input:
        "{sample}.fasta",
    output:
        "{sample}.fasta.fai",
    wrapper:
        "v1.21.2/bio/samtools/faidx"


"""
rule get_single_chromosome_reference:
    input:
        ref="data/{build}/reference.fa",
        index="data/{build}/reference.fa.fai",
    output:
        "data/{build}/raw/{build}_{chromosome}.fa"
    params:
        "{chromosome}"
    wrapper:
        "v1.21.2/bio/samtools/faidx"
"""

rule get_dataset_reference:
    input:
        ref = "data/{genome_build}/reference.fa",
        index = "data/{genome_build}/reference.fa.fai",
    output:
        ReferenceGenome.path()
        #"data/{reference}/{individual}/{size}/reference.fa"
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