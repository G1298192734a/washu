#!/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3
from pybedtools import BedTool as bt
import numpy as np

bed=bt("/work/frasergen/3D/work/shaojie/script/HiC/WGBS/methyl_snake/01.ref/cpg.bed")
med_len=np.mean([i.length for i in bed])
body_len=(med_len//20+1)*20
bin_size=20#med_len//20+1




