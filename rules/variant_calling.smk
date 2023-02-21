# Variant calling is done on a single chromosome to save time
# but is done on reads mapped to/from the whole genome to not introduce any bias

def get_variant_calling_chromosome(wildcards):
    chromosomes = config["genomes"][wildcards.genome_build][wildcards.individual][wildcards.dataset_size]["chromosomes"].split(",")
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
    threads: 2
    wrapper:
        "v1.21.4/bio/sambamba/sort"


rule filter_bam_on_mapq:
    input:
        f"data/{parameters.until('n_threads')}/mapped.bam"
    output:
        f"data/{parameters.until('min_mapq')}/mapped.bam"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools view -b -q {wildcards.min_mapq} {input} > {output} 
        """


rule index_bam:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools index {input}"


rule call_variants:
    input:
        sorted_bam = f"data/{parameters.until('min_mapq')}/mapped.sorted.bam",
        bamindex = f"data/{parameters.until('min_mapq')}/mapped.sorted.bam.bai",
        reference = f"data/{parameters.until('dataset_size')}/reference.fa",
        reference_index = f"data/{parameters.until('dataset_size')}/reference.fa.fai",
        reference_dict = f"data/{parameters.until('dataset_size')}/reference.dict",
    output:
        f"data/{parameters.until('min_mapq')}/variant_calls.vcf.gz"
    params:
        chromosome=get_variant_calling_chromosome
    threads:
        1
    conda:
        "../envs/gatk.yml"
    shell:
        "gatk HaplotypeCaller "
        "--reference {input.reference} "
        "--input {input.sorted_bam} "
        "--native-pair-hmm-threads 1 "
        "--output {output} "
        "--intervals {params.chromosome} "


rule make_truth_vcf_for_chromosome:
    input:
        variants="data/{genome_build}/{individual}/variants.vcf.gz",
        index="data/{genome_build}/{individual}/variants.vcf.gz.tbi",
    output:
        "data/{genome_build,\w+}/{individual}/{dataset_size}/variant_calling_truth.vcf.gz"
    params:
        chromosome=get_variant_calling_chromosome
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools view --regions {params.chromosome} -O z {input.variants} > {output}
        """



rule run_happy:
    input:
        variant_calls=f"data/{parameters.until('min_mapq')}/variant_calls.vcf.gz",
        truth_vcf=f"data/{parameters.until('dataset_size')}/variant_calling_truth.vcf.gz",
        truth_regions=f"data/{parameters.until('individual')}/truth_regions.bed",
        ref=f"data/{parameters.until('genome_build')}/reference.fa"
    output:
        f"data/{parameters.until('min_mapq')}/happy.extended.csv",
        f"data/{parameters.until('min_mapq')}/happy.summary.csv",
    params:
        output_base_name=lambda wildcards, input, output: output[0].split(".")[0]
    conda:
        "../envs/happy.yml"
    shell:
        """
        export HGREF={input.ref} && \
        hap.py {input.truth_vcf} {input.variant_calls} \
        --no-leftshift \
        -r {input.ref} \
        -o {params.output_base_name} \
        -f {input.truth_regions} \
        --no-decompose --engine=vcfeval 
        """


"""
rule run_happy2:
    input:
        truth="data/{genome_build}/{individual}/{size}/variant_calling_truth.vcf.gz",
        query="data/{genome_build}/{individual}/{size}/{config}/variant_calls.vcf.gz",
        truth_regions="data/{genome_build}/{individual}/truth_regions.bed",
        genome="data/{genome_build}/reference.fa",
        genome_index="data/{genome_build}/reference.fa.fai"
    output:
        multiext("data/{genome_build,\w+}/{individual,\w+}/{size,\w+}/{config}/happy",".runinfo.json",".vcf.gz",".summary.csv",
            ".extended.csv",".metrics.json.gz",".roc.all.csv.gz",
            ".roc.Locations.INDEL.csv.gz",".roc.Locations.INDEL.PASS.csv.gz",
            ".roc.Locations.SNP.csv.gz",".roc.tsv")
    params:
        engine="vcfeval",
        prefix=lambda wc, input, output: output[0].split('.')[0],
        ## parameters such as -L to left-align variants
        extra="--verbose"
    threads: 2
    wrapper: "v1.23.4/bio/hap.py/hap.py"
"""

rule get_variant_calling_result:
    input:
        f"data/{parameters.until('min_mapq')}/happy.summary.csv"
    output:
        f"data/{parameters}/variant_calling_{{type, recall|one_minus_precision|f1score}}.txt"
    run:
        import pandas as pd
        data = pd.read_csv(input[0])

        index = 2
        if wildcards.variant_calling_type == "indels":
            index = 0

        names = {"recall": "METRIC.Recall",
                 "one_minus_precision": "METRIC.Precision",
                 "f1score": "METRIC.F1_Score"}

        result = data.iloc[index][names[wildcards.type]]
        with open(output[0], "w") as f:
            if wildcards.type == "one_mnus_precision":
                result = 1 - result

            f.write(str(result) + "\n")

