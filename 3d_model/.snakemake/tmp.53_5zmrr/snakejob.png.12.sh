#!/bin/sh
# properties = {"type": "single", "rule": "png", "local": false, "input": ["/work/frasergen/backup/3d/project/Interactome/240033_kcybgj_20240313/04.HiC/14.3d_model/02.3d_model/Delta-abrB/Delta-abrB.pdb", "/work/frasergen/backup/3d/project/Interactome/240033_kcybgj_20240313/04.HiC/14.3d_model/02.3d_model/Delta-abrB/Delta-abrB.pml"], "output": ["/work/frasergen/backup/3d/project/Interactome/240033_kcybgj_20240313/04.HiC/14.3d_model/03.gif/Delta-abrB/pngs/fig_090.png"], "wildcards": {"sample": "Delta-abrB", "batch": "090"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 12, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/3d_model && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3 \
-m snakemake /work/frasergen/backup/3d/project/Interactome/240033_kcybgj_20240313/04.HiC/14.3d_model/03.gif/Delta-abrB/pngs/fig_090.png --snakefile /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.53_5zmrr /work/frasergen/backup/3d/project/Interactome/240033_kcybgj_20240313/04.HiC/14.3d_model/02.3d_model/Delta-abrB/Delta-abrB.pdb /work/frasergen/backup/3d/project/Interactome/240033_kcybgj_20240313/04.HiC/14.3d_model/02.3d_model/Delta-abrB/Delta-abrB.pml --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid.yaml -p --allowed-rules png --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.53_5zmrr/12.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.53_5zmrr/12.jobfailed; exit 1)

