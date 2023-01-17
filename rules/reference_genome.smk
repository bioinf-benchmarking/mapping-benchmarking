

rule download_reference:
    output:
        "data/reference_genomes/raw/{reference}.2bit"
    shell:
        "wget -O {output} https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.reference}/bigZips/{wildcards.reference}.2bit"


rule convert_reference_genome_to_fasta:
    input:
        "data/reference_genomes/raw/{sample}.2bit"
    output:
        "data/reference_genomes/raw/{sample}.fa"
    log:
        "logs/{sample}.2bit_to_fa.log"
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


rule get_single_chromosome_reference:
    input:
        ref="data/reference_genome/raw/{build}.fa",
        index="data/reference_genome/raw/{build}.fa.fai",
    output:
        "data/reference_genome/raw/{build}_{chromosome}.fa"
    params:
        "{chromosome}"
    wrapper:
        "v1.21.2/bio/samtools/faidx"


rule get_dataset_reference:
    input:
        ref=lambda wildcards: "data/reference_genomes/raw/" + config["simulations"][wildcards.dataset]["genome"]  + ".fa",
        index=lambda wildcards: "data/reference_genomes/raw/" + config["simulations"][wildcards.dataset]["genome"]  + ".fa.fai",
    output:
        "data/reference_genomes/{dataset}.fa"
    conda:
        "envs/samtools.yml"
    params:
        regions=lambda wildcards: config["simulations"][wildcards.dataset]["chromosomes"].replace(",", " ")
    shell:
        "samtools faidx {input.ref} {params.regions} > {output}"
    #wrapper:
    #    "v1.21.2/bio/samtools/faidx"
