
def get_truth_vcf_command(wildcards, input, output):
    individual_config = config["genomes"][wildcards.genome_build][wildcards.individual]
    if individual_config["simulated"]:
        tmp_output = output.vcf.split(".gz")[0]
        parameters = ' '.join(config["mason_variator_parameters"].split("\n"))
        return f"mason_variator --seed 123 -ir {input.reference} -ov {tmp_output} {parameters} && cat {tmp_output} | python scripts/assign_random_genotypes_to_vcf.py | bgzip -c > {output.vcf}"
    else:
        url = individual_config["vcf_url"]
        return f"wget -O {output} {url}"


rule get_truth_vcf:
    input:
        reference=GenomeBuild.path() + "/reference.fa",
        fai=GenomeBuild.path() + "/reference.fa.fai",
    output:
        vcf=Individual.path() + "/variants.vcf.gz"
    params:
        command=get_truth_vcf_command
    conda:
        "../envs/mason.yml"
    shell:
        "{params.command}"
