# Assigns numeric incremental IDs to a fq file

import sys
read_id = 0
for i, line in enumerate(sys.stdin):
    if i % 4 == 0:
        print(f"@{read_id}")
        read_id += 1
    else:
        print(line.strip())