#!/bin/sh
# properties = {"type": "single", "rule": "single", "local": false, "input": ["/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/02.single_pre/CapHiC-test", "/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/02.single_pre/CapHiC-test/work.sh"], "output": ["/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/CapHiC-test", "/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/CapHiC-test/work.sh"], "wildcards": {"spl": "CapHiC-test"}, "params": {}, "log": ["logs/single_CapHiC-test.log"], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 4, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/washU/upload && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3 \
-m snakemake /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/CapHiC-test --snakefile /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/upload_hic3.0_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.5fal_0wm /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/02.single_pre/CapHiC-test /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/analysis20240624_washu/02.single_pre/CapHiC-test/work.sh --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/upload_hic3.yaml -p --allowed-rules single --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.5fal_0wm/4.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.5fal_0wm/4.jobfailed; exit 1)

