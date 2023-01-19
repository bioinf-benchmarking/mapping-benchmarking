# Assigns numeric incremental IDs to a fq file

import sys
read_id = 0
for i, line in enumerate(sys.stdin):
    if line.startswith("@"):
        continue

    line = line.strip().split("\t")
    line[0] = str(read_id)
    print("\t".join(line))
    read_id += 1