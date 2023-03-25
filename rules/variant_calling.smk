from mapping_benchmarking.config import WholeGenomeMappedReads, FilteredWholeGenomeMappedReads, VariantCalls, \
    ReferenceGenome, MapQFilteredWholeGenomeMappedReads, Individual, VariantCallingAccuracy


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
        "{some_file}.bam"
    output:
        "{some_file}.sorted.bam"
    params:
        ""  # optional parameters
    threads: 2
    wrapper:
        "v1.21.4/bio/sambamba/sort"


rule filter_bam_on_mapq:
    input:
        WholeGenomeMappedReads.path()
    output:
        MapQFilteredWholeGenomeMappedReads.path()
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
        sorted_bam = MapQFilteredWholeGenomeMappedReads.path(file_ending = ".sorted.bam"),
        bamindex = MapQFilteredWholeGenomeMappedReads.path(file_ending = ".sorted.bam.bai"),
        reference = ReferenceGenome.path(),
        reference_index = ReferenceGenome.path(file_ending=".fa.fai"),
        reference_dict= ReferenceGenome.path(file_ending=".dict"),
    output:
        VariantCalls.path()
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
        variants=Individual.path() + "/variants.vcf.gz",
        index=Individual.path() + "/variants.vcf.gz.tbi"
    output:
        ReferenceGenome.path(file_ending="") + "/variant_calling_truth.vcf.gz"
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
        variant_calls=VariantCalls.path(),
        truth_vcf=ReferenceGenome.path(file_ending="") + "/variant_calling_truth.vcf.gz",
        truth_regions= Individual.path() + "/truth_regions.bed",
        ref=ReferenceGenome.path()
    output:
        VariantCalls.path(file_ending="") + "/happy.extended.csv",
        VariantCalls.path(file_ending="") + "/happy.summary.csv",
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


rule get_variant_calling_result:
    input:
        VariantCalls.path(file_ending="") + "/happy.summary.csv",
    output:
        VariantCallingAccuracy.path()
    run:
        import pandas as pd
        data = pd.read_csv(input[0])

        index = 2
        if wildcards.variant_calling_type == "indels":
            index = 0

        names = {"VariantCallingRecall": "METRIC.Recall",
                 "VariantCallingOneMinusPreicision": "METRIC.Precision",
                 "VariantCallingF1Score": "METRIC.F1_Score"}

        result = data.iloc[index][names[wildcards.type]]
        with open(output[0], "w") as f:
            if wildcards.type == "VariantCallingOneMinusPrecision":
                result = 1 - result

            f.write(str(result) + "\n")
