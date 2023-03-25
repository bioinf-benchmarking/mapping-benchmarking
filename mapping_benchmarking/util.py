from .config import GenericReads

def get_input_reads(wildcards):
    file = "reads.fq.gz" if not "paired_end" in wildcards.read_config else ["reads1.fq.gz", "reads2.fq.gz"]
    return GenericReads.path(file=file)
