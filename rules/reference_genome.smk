

rule download_reference:
    output:
        "data/{reference}/reference.2bit"
    shell:
        "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.reference}/bigZips/{wildcards.reference}.2bit"


rule convert_reference_genome_to_fasta:
    input:
        "data/{sample}/reference.2bit"
    output:
        "data/{sample,\w+}/reference.fa"
    wrapper:
        "v1.21.2/bio/ucsc/twoBitToFa"


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
        ref = "data/{reference}/reference.fa",
        index = "data/{reference}/reference.fa.fai",
        #ref=lambda wildcards: "data/reference_genomes/raw/" + config["variant_sources"][config["individuals"][wildcards.individual]["variant_source"]]["genome"]  + ".fa",
        #index=lambda wildcards: "data/reference_genomes/raw/" + config["variant_sources"][config["individuals"][wildcards.individual]["variant_source"]]["genome"]  + ".fa.fai",
    output:
        "data/{reference}/{individual}/{size}/reference.fa"
    conda:
        "../envs/samtools.yml"
    params:
        regions=lambda wildcards: config["genomes"][wildcards.reference][wildcards.individual][wildcards.size]["chromosomes"].replace(",", " ")
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