#!/bin/sh
# properties = {"type": "single", "rule": "matrix", "local": false, "input": [], "output": ["/work/frasergen/backup/3d/project/Interactome/232385_Schizonepeta_tenuifolia/05.personal/analysis_3D_model_20240507/01.matrix_fill/Root/Root_filled_matrix.txt"], "wildcards": {"sample": "Root"}, "params": {"hm": "/work/frasergen/backup/3d/project/Interactome/232385_Schizonepeta_tenuifolia/03.HiC/04.matrix/Root.100000.cool", "python": "/public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin/python", "hm2matrix": "/public/frasergen/3D/work/wanghao/my_scripts/python/202212/hiclib_hm2txt/step01.txt.py", "cool2matrix": "/work/frasergen/backup/3d/project/Interactome/221514-01_zhenjun_20240103/07.personal/analysis_20240205_3d_model/genome_cool_to_dense.py", "matrix_fill": "/work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_soft/fill_blanks_v0.1.awk"}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 5, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/3d_model && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3.7 \
-m snakemake /work/frasergen/backup/3d/project/Interactome/232385_Schizonepeta_tenuifolia/05.personal/analysis_3D_model_20240507/01.matrix_fill/Root/Root_filled_matrix.txt --snakefile /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_snake_test.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.h86y3725 --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid1.yaml -p --allowed-rules matrix --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.h86y3725/5.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/3d_model/.snakemake/tmp.h86y3725/5.jobfailed; exit 1)

