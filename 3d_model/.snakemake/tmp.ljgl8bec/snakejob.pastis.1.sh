#!/bin/sh
# properties = {"type": "single", "rule": "pastis", "local": false, "input": ["/work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/01.matrix_fill/AA"], "output": ["/work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/02.3d_model/AA/AA.pdb"], "wildcards": {"sample": "AA"}, "params": {"python": "/public/frasergen/3D/work/wanghao/myconda/pastis-0.1.0/bin/python", "heat": "/public/frasergen/3D/work/wanghao/pipeline/3D_model_pastis-0.1.0/3d.py", "dir": "/work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/02.3d_model/AA", "pdb": "/work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/02.3d_model/AA/MDS.test.fa.fai-100000.pdb", "filled": "/work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/01.matrix_fill/AA/AA_corrected_matrix.txt.filled_0.90"}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 1, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/3d_model && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3.7 \
-m snakemake /work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/02.3d_model/AA/AA.pdb --snakefile /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_test_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.ljgl8bec /work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/01.matrix_fill/AA --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid.yaml -p --allowed-rules pastis --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.ljgl8bec/1.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.ljgl8bec/1.jobfailed; exit 1)

