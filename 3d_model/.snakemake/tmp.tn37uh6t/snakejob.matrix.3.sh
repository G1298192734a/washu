#!/bin/sh
# properties = {"type": "single", "rule": "matrix", "local": false, "input": [], "output": ["/work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis_chr20_3d_model_20240531/01.matrix_fill/MH-zm/MH-zm_filled_matrix.txt"], "wildcards": {"sample": "MH-zm"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 3, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/3d_model && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3 \
-m snakemake /work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis_chr20_3d_model_20240531/01.matrix_fill/MH-zm/MH-zm_filled_matrix.txt --snakefile /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.tn37uh6t --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid.yaml -p --allowed-rules matrix --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.tn37uh6t/3.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.tn37uh6t/3.jobfailed; exit 1)

