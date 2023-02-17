

rule strobealign_map:
  input:
    reads=get_input_reads,
    ref=f"data/{parameters.until('dataset_size')}/reference.fa"
  output:
      f"data/{parameters.until('n_threads')(method='strobealign')}/mapped.bam"
  conda:
      "../envs/strobealign.yml"
  threads: lambda wildcards: int(wildcards.n_threads)
  benchmark:
      f"data/{parameters.until('n_threads')(method='strobealign')}/benchmark.csv"
  shell:
      "mkdir -p dirname {output} && "
      "strobealign -t {wildcards.n_threads} {input.ref} {input.reads} | "
      "samtools view -b -h - | "
      # add read group (needed for variant calling)
      "samtools addreplacerg -r 'ID:sample\\tSM:sample' -O bam -o {output} - "


