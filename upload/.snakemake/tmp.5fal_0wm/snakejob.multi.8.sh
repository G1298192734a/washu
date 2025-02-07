#!/bin/sh
# properties = {"type": "single", "rule": "multi", "local": false, "input": ["/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/03.multi_pre/CapHiC-test.minus.CapHiC-control", "/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/03.multi_pre/CapHiC-test.minus.CapHiC-control/work.sh"], "output": ["/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/CapHiC-test.minus.CapHiC-control", "/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/CapHiC-test.minus.CapHiC-control/work.sh"], "wildcards": {"minus": "CapHiC-test.minus.CapHiC-control"}, "params": {}, "log": ["logs/multi_CapHiC-test.minus.CapHiC-control.log"], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 8, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/washU/upload && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3 \
-m snakemake /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/CapHiC-test.minus.CapHiC-control --snakefile /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/upload_hic3.0_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.5fal_0wm /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/03.multi_pre/CapHiC-test.minus.CapHiC-control /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/03.multi_pre/CapHiC-test.minus.CapHiC-control/work.sh --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/upload_hic3.yaml -p --allowed-rules multi --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.5fal_0wm/8.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.5fal_0wm/8.jobfailed; exit 1)

