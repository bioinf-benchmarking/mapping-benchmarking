# Assigns numeric incremental IDs to a fq file
import logging
logging.basicConfig(level=logging.INFO)
import sys

read_id = 0
is_paired_end = False
i = 0
for line in sys.stdin:
    if line.startswith("@"):
        continue

    line = line.strip().split("\t")
    if line[6] != "*":
        if not is_paired_end:
            logging.info("Assuming reads are paired end")
        is_paired_end = True

    read_name = str(read_id)
    if is_paired_end:
        read_name += "/" + str((i % 2)+1)

    line[0] = str(read_name)
    print("\t".join(line))

    if not is_paired_end or i % 2 == 1:
        read_id += 1

    i += 1


logging.info("In total %d reads" % i)