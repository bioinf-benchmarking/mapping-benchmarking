# Assigns numeric incremental IDs to a fq file
import logging
logging.basicConfig(level=logging.INFO)
import sys

is_paired_end = False
read_id = 0
for i, line in enumerate(sys.stdin):
    if i % 4 == 0:
        if "/" in line:
            is_paired_end = True

        name = str(read_id)
        if is_paired_end:
            name += "/" + str((i % 8) // 4 + 1)

        print(f"@{name}")

        if not is_paired_end or i % 8 == 4:
            read_id += 1
    else:
        print(line.strip())