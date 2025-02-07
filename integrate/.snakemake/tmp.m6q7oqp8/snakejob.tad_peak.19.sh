#!/bin/sh
# properties = {"type": "single", "rule": "tad_peak", "local": false, "input": [], "output": [], "wildcards": {"diff": "CapHiC-test.minus.CapHiC-control", "atac": "HiC_CUT_H3K27me3"}, "params": {}, "log": ["/work/frasergen/3D/work/shaojie/script/HiC/integrate/analysis_HiC_ATAC_RNA/01.integrate/CapHiC-test.minus.CapHiC-control/HiC_CUT_H3K27me3/04.tad_peak/work.sh"], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 19, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/integrate && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3.7 \
-m snakemake /work/frasergen/3D/work/shaojie/script/HiC/integrate/analysis_HiC_ATAC_RNA/01.integrate/CapHiC-test.minus.CapHiC-control/HiC_CUT_H3K27me3/04.tad_peak/work.sh --snakefile /work/frasergen/3D/work/shaojie/script/HiC/integrate/integrate_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.m6q7oqp8 --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/integrate/integrate.yaml -p --allowed-rules tad_peak --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.m6q7oqp8/19.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.m6q7oqp8/19.jobfailed; exit 1)

