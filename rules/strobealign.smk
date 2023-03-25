from mapping_benchmarking.config import *


# seperate rule so that indexing time is not included in benchmark
# index depends on reads, but output file name is not deterministic, so create a checkpoint file that
# depends on the reads
rule strobealign_index:
    input:
        ref=ReferenceGenome.path(file_ending=".fa"),
        reads=get_input_reads
    output:
        touch(ReferenceGenome.path(file_ending="") + "/{read_config}/strobealign-index-created")
        #touch(f"data/{reference_genome}/{{config}}/strobealign-index-created")
    conda:
        "../envs/strobealign.yml"
    shell:
        "strobealign --create-index {input.ref} {input.reads}"


rule strobealign_map:
  input:
    reads = get_input_reads,
    ref = ReferenceGenome.path(file_ending=".fa"),
    #ref=f"data/{reference_genome}/reference.fa",
    checkpoint_index_is_created = ReferenceGenome.path(file_ending="") + "/{read_config}/strobealign-index-created"
  output:
      #f"data/{reference_genome}/{{config}}/strobealign/{{n_threads}}/without_readgroup.bam"
      reads = GenericMappedReads.as_output(method='strobealign')
  conda:
      "../envs/strobealign.yml"
  threads: lambda wildcards: int(wildcards.n_threads)
  benchmark:
      GenericMappedReads.as_output(method='strobealign',file_ending=".benchmark.csv")
  shell:
      "mkdir -p $(dirname {output}) && "
      "strobealign --rg SM:sample --rg-id sample -t {wildcards.n_threads} --use-index {input.ref} {input.reads} | "
      "samtools view -o {output} -"
