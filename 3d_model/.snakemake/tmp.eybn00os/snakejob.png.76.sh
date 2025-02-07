#!/bin/sh
# properties = {"type": "single", "rule": "png", "local": false, "input": ["/work/frasergen/backup/3d/project/Interactome/240137_YanDanBaoJun_20240425/06.personal/analysis_3D_model_20240524/02.3d_model/TD_60MM12h/TD_60MM12h.pdb", "/work/frasergen/backup/3d/project/Interactome/240137_YanDanBaoJun_20240425/06.personal/analysis_3D_model_20240524/02.3d_model/TD_60MM12h/TD_60MM12h.pml"], "output": ["/work/frasergen/backup/3d/project/Interactome/240137_YanDanBaoJun_20240425/06.personal/analysis_3D_model_20240524/03.gif/TD_60MM12h/pngs/fig_150.png"], "wildcards": {"sample": "TD_60MM12h", "batch": "150"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 76, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/3d_model && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3 \
-m snakemake /work/frasergen/backup/3d/project/Interactome/240137_YanDanBaoJun_20240425/06.personal/analysis_3D_model_20240524/03.gif/TD_60MM12h/pngs/fig_150.png --snakefile /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.eybn00os /work/frasergen/backup/3d/project/Interactome/240137_YanDanBaoJun_20240425/06.personal/analysis_3D_model_20240524/02.3d_model/TD_60MM12h/TD_60MM12h.pdb /work/frasergen/backup/3d/project/Interactome/240137_YanDanBaoJun_20240425/06.personal/analysis_3D_model_20240524/02.3d_model/TD_60MM12h/TD_60MM12h.pml --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid.yaml -p --allowed-rules png --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.eybn00os/76.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.eybn00os/76.jobfailed; exit 1)

