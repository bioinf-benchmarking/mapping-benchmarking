

rule minimap_index:
    input:
        "{data}/{ref}.fa",
    output:
        "{data}/{ref}.mmi"
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -ax sr -d {output} {input}
        """

rule minimap_map:
    input:
        target = ReferenceGenome.path(file_ending=".mmi"),
        #target = f"data/{reference_genome}/reference.mmi",
        query = get_input_reads  # "{data}/whole_genome_single_end/{config}/reads.fq.gz",
    output:
        #"{data}/whole_genome_single_end/{config}/minimap/{n_threads}/mapped.bam",
        reads=GenericMappedReads.as_output(method='minimap')
        #f"data/{reference_genome}/{{config}}/minimap/{{n_threads}}/mapped.bam"
    params:
        extra = r"-ax sr -R '@RG\tID:sample\tSM:sample' -a",  # optional
        sorting = "none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra = "",  # optional: extra arguments for samtools/picard

    benchmark:
        #"{data}/whole_genome_single_end/{config}/minimap/{n_threads}/benchmark.csv",
        #f"data/{reference_genome}/{{config}}/minimap/{{n_threads}}/benchmark.csv"
        GenericMappedReads.as_output(method='minimap', file_ending=".benchmark.csv")
    threads: lambda wildcards: int(wildcards.n_threads)
    conda:
        "../envs/minimap2.yml"
    shell:
        # hack with removing column 16, gatk does not handle the tp:A:P field
        """
        minimap2 -t {wildcards.n_threads} --MD -ax sr -R '@RG\\tID:sample\\tSM:sample' -a {input} | cut --complement -f 17 | samtools view -b -h - > {output}
        """
