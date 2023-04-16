# mason_variator outputs vcf without genotypes
# this script just assign a random genotype to each variant

import sys
import random

genotypes = ["0|0", "0|1", "1|0", "1|1"]

for line in sys.stdin:
    if line.startswith("#"):
        print(line.strip())
        continue

    l = line.split()
    l[8] = "GT"
    l[9] = random.choice(genotypes)
    print("\t".join(l).strip())