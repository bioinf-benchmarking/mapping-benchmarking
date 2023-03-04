from shared_memory_wrapper import from_file

coordinate_map = from_file(snakemake.input.coordinate_map)
out_file = open(snakemake.output.sam, "w")
out_file_txt = open(snakemake.output.txt, "w")

for line in open(snakemake.input.truth_positions):
    if line.startswith("@"):
        out_file.write(line)
        continue

    l = line.split()
    chromosome = l[2]
    start = int(l[3])
    length = len(l[9])
    end = start + length
    n_variants = 0
    if coordinate_map.haplotype_has_variant_between(chromosome, start, end):
        n_variants = 1

    #line = line.strip() + "\tNVARIANTS:i:" + str(n_variants) + "\n"
    out_file.write(line)
    out_file_txt.write(str(n_variants) + "\n")
