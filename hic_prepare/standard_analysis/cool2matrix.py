#!/public/frasergen/PUB/software/python/Python-3.7.7/bin/python3
import cooler,sys
import numpy as np
from diff_matrix import parse_genome_cool

if __name__ == "__main__":
  parse_genome_cool(sys.argv[1], sys.argv[2])


