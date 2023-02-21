
# seperate rule so that indexing time is not included in benchmark
# index depends on reads, but output file name is not deterministic, so create a checkpoint file that
# depends on the reads
rule strobealign_index:
    input:
        ref=f"data/{parameters.until('dataset_size')}/reference.fa",
        reads=get_input_reads
    output:
        touch(f"data/{parameters.until('n_reads')}/strobealign-index-created")
    conda:
        "../envs/strobealign.yml"
    shell:
        "strobealign --create-index {input.ref} {input.reads}"


rule strobealign_map:
  input:
    reads=get_input_reads,
    ref=f"data/{parameters.until('dataset_size')}/reference.fa",
    checkpoint_index_is_created=f"data/{parameters.until('n_reads')}/strobealign-index-created"
  output:
      f"data/{parameters.until('n_threads')(method='strobealign')}/without_readgroup.bam"
  conda:
      "../envs/strobealign.yml"
  threads: lambda wildcards: int(wildcards.n_threads)
  benchmark:
      f"data/{parameters.until('n_threads')(method='strobealign')}/benchmark.csv"
  shell:
      "mkdir -p $(dirname {output}) && "
      "strobealign -t {wildcards.n_threads} --use-index {input.ref} {input.reads} | "
      "samtools view -o {output} -"



# have this as separate rule to not be includeded in benchmark time
rule add_read_group_to_strobealign:
    input:
        f"data/{parameters.until('n_threads')(method='strobealign')}/without_readgroup.bam"
    output:
        f"data/{parameters.until('n_threads')(method='strobealign')}/mapped.bam"
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools addreplacerg -r 'ID:sample\\tSM:sample' -O bam -o {output} {input}"