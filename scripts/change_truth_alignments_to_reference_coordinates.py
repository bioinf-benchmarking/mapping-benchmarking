from shared_memory_wrapper import from_file

coordinate_map = from_file(snakemake.input.coordinate_map)
out_file = open(snakemake.output[0], "w")
n = 0
for line in open(snakemake.input.truth_positions):
    if line.startswith("@"):
        out_file.write(line)
    else:
        l = line.split("\t")
        chromosome = l[2]
        position = int(l[3])
        ref_position = coordinate_map.convert(chromosome, position)
        l[3] = str(ref_position)
        out_file.write("\t".join(l))
        n += 1

print("%d reads processed" % n)
