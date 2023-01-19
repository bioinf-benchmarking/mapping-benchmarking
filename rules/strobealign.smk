

rule strobealign_map:
  input:
    reads="data/simulated_reads/{dataset}/simulated_reads.fq.gz",
    ref="data/reference_genomes/{dataset}.fa"
  output:
      sam="data/mapping/strobealign/{dataset}.sam",
  conda:
      "../envs/strobealign.yml"
  shell:
      "strobealign {input.ref} {input.reads} > {output.sam}"


