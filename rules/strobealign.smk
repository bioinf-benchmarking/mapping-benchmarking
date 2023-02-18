
# seperate rule so that indexing time is not included in benchmark
# index depends on read length and thus the reads
rule strobealign_index:
    input:
        ref=f"data/{parameters.until('dataset_size')}/reference.fa"
    output:
        f"data/{parameters.until('dataset_size')}/reference.fa.r{{r}}.sti"
    conda:
        "../envs/strobealign.yml"
    shell:
        "strobealign --create-index -r {wildcards.r} {input.ref}"


def strobealign_r(read_length):
    return 50 * round(int(read_length) / 50)


rule strobealign_map:
  input:
    reads=get_input_reads,
    ref=f"data/{parameters.until('dataset_size')}/reference.fa",
    index=lambda wildcards: f"data/{parameters.until('dataset_size')}/reference.fa.r{strobealign_r(wildcards.read_length)}.sti"
  output:
      f"data/{parameters.until('n_threads')(method='strobealign')}/without_readgroup.bam"
  conda:
      "../envs/strobealign.yml"
  threads: lambda wildcards: int(wildcards.n_threads)
  params:
      # strobealign rounds r to nearest 50, so we need to do the same to know what the index neede will be
      r=lambda wildcards: strobealign_r(wildcards.read_length)
  benchmark:
      f"data/{parameters.until('n_threads')(method='strobealign')}/benchmark.csv"
  shell:
      "mkdir -p $(dirname {output}) && "
      "strobealign -t {wildcards.n_threads} -r {params.r} --use-index {input.ref} {input.reads} | "
      "samtools view -b -h - > {output} "



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
