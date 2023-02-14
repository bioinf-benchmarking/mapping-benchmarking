

rule minimap_index:
    input:
        "{data}/reference.fa",
    output:
        "{data}/reference.mmi"
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -ax sr -d {output} {input}
        """

rule minimap_single_end:
    input:
        target = "{data}/reference.mmi",
        query = "{data}/whole_genome_single_end/{config}/reads.fq.gz",
    output:
        "{data}/whole_genome_single_end/{config}/minimap/{n_threads}/mapped.bam",
    params:
        extra = r"-ax sr -R '@RG\tID:sample\tSM:sample' -a",  # optional
        sorting = "none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra = "",  # optional: extra arguments for samtools/picard

    benchmark:
        "{data}/whole_genome_single_end/{config}/minimap/{n_threads}/benchmark.csv",
    threads: lambda wildcards: int(wildcards.n_threads)
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -t {wildcards.n_threads} --MD -ax sr -R '@RG\\tID:sample\\tSM:sample' -a {input} | samtools view -b -h - > {output}
        """
    #wrapper:
    #    "v1.21.2/bio/minimap2/aligner"


rule minimap_paired_end:
    input:
        target = "{data}/reference.mmi",
        query = ["{data}/whole_genome_paired_end/{config}/reads" + n + ".fq.gz" for n in ("1", "2")],
    output:
        "{data}/whole_genome_paired_end/{config}/minimap/{n_threads}/mapped.bam",
    params:
        extra=r"-ax sr -R '@RG\tID:sample\tSM:sample' -a",# optional
        sorting="none",# optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",# optional: extra arguments for samtools/picard
    benchmark:
        "{data}/whole_genome_paired_end/{config}/minimap/{n_threads}/benchmark.csv",
    threads: lambda wildcards: int(wildcards.n_threads)
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -t {wildcards.n_threads} --MD -ax sr -R '@RG\\tID:sample\\tSM:sample' -a {input} | samtools view -b -h - > {output}
        """
    #wrapper:
    #    "v1.21.2/bio/minimap2/aligner"
