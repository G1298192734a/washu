#!/bin/sh
# properties = {"type": "single", "rule": "png", "local": false, "input": ["/work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/03.gif/AC/pmls/gif_009.pml"], "output": ["/work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/03.gif/AC/pngs/fig_009.png"], "wildcards": {"sample": "AC", "batch": "009"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 139, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/3d_model && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3.7 \
-m snakemake /work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/03.gif/AC/pngs/fig_009.png --snakefile /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_test_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.h2giwvet /work/frasergen/3D/work/shaojie/script/HiC/3d_model/analysis20240419_3D_model/03.gif/AC/pmls/gif_009.pml --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid.yaml -p --allowed-rules png --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.h2giwvet/139.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.h2giwvet/139.jobfailed; exit 1)

