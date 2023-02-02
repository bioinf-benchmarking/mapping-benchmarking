# Variant calling is done on a single chromosome to save time
# but is done on reads mapped to/from the whole genome to not introduce any bias

def get_variant_calling_chromosome(wildcards):
    chromosomes = config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.size]["chromosomes"].split(",")
    return chromosomes[0]


rule make_dict_file:
    input: "{file}.fa"
    output: "{file}.dict"
    #conda: "../envs/picard.yml"
    wrapper:
        "v1.21.4/bio/picard/createsequencedictionary"

rule sambamba_sort:
    input:
        "{dir}/mapped.bam"
    output:
        "{dir}/mapped.sorted.bam"
    params:
        ""  # optional parameters
    threads: 8
    wrapper:
        "v1.21.4/bio/sambamba/sort"


rule call_variants:
    input:
        sorted_bam = "data/{genome_build}/{individual}/{size}/whole_genome_{pairing}_end/{config}/mapped.sorted.bam",
        reference = "data/{genome_build}/{individual}/{size}/reference.fa",
        reference_index = "data/{genome_build}/{individual}/{size}/reference.fa.fai",
        reference_dict = "data/{genome_build}/{individual}/{size}/reference.dict",
    output:
        "{data}/{genome_build}/{individual}/{size}/whole_genome_{pairing}_end/{config}/variant_calls.vcf.gz"
    params:
        chromosome=get_variant_calling_chromosome
    conda:
        "../envs/gatk.yml"
    shell:
        "gatk HaplotypeCaller "
        "--reference {input.reference} "
        "--input {input.sorted_bam} "
        "--native-pair-hmm-threads 1 "
        "--output {output} "
        "--intervals {params.chromosome} "
        "--minimum-mapping-quality 20 "


rule make_truth_vcf_for_chromosome:
    input:
        variants="data/{genome_build}/{individual}/variants.vcf.gz",
        index="data/{genome_build}/{individual}/variants.vcf.gz.tbi",
    output:
        "data/{genome_build,\w+}/{individual}/{size}/variant_calling_truth.vcf.gz"
    params:
        chromosome=get_variant_calling_chromosome
    shell:
        """
        bcftools view --regions {params.chromosome} -O z {input.variants} > {output}
        """


rule run_happy:
    input:
        variant_calls="data/{genome_build}/{individual}/{size}/{config}/variant_calls.vcf.gz",
        truth_vcf="data/{genome_build}/{individual}/{size}/variant_calling_truth.vcf.gz",
        truth_regions="data/{genome_build}/{individual}/truth_regions.bed",
        ref="data/{genome_build}/reference.fa"
    output:
        "data/{genome_build,\w+}/{individual,\w+}/{size,\w+}/{config}/happy.extended.csv"
    conda:
        "../envs/happy.yml"
    shell:
        """
        hap.py {input.truth_vcf} {input.variant_calls} \
        --no-leftshift \
        -r {input.ref} \
        -o data/{wildcards.genome_build}/{wildcards.individual}/{wildcards.size}/{wildcards.config}/happy \
        -f {input.truth_regions} \
        --no-decompose --engine=vcfeval 
        """

