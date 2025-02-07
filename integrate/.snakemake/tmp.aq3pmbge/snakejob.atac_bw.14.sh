#!/bin/sh
# properties = {"type": "single", "rule": "atac_bw", "local": false, "input": ["/work/frasergen/backup/3d/project/CUT/2024/231661_ren/CTCF/01.QC/01.analysis/CTCF-test/CTCF-test.CPM.bw", "/work/frasergen/backup/3d/project/CUT/2024/231661_ren/CTCF/01.QC/01.analysis/CTCF-control/CTCF-control.CPM.bw"], "output": ["/work/frasergen/3D/work/shaojie/script/HiC/integrate/analysis_HiC_ATAC_RNA/00.pre/HiC_CUT_CTCF_pre/02.cpm_bw"], "wildcards": {"atac": "HiC_CUT_CTCF"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 14, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/integrate && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3.7 \
-m snakemake /work/frasergen/3D/work/shaojie/script/HiC/integrate/analysis_HiC_ATAC_RNA/00.pre/HiC_CUT_CTCF_pre/02.cpm_bw --snakefile /work/frasergen/3D/work/shaojie/script/HiC/integrate/integrate_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.aq3pmbge /work/frasergen/backup/3d/project/CUT/2024/231661_ren/CTCF/01.QC/01.analysis/CTCF-test/CTCF-test.CPM.bw /work/frasergen/backup/3d/project/CUT/2024/231661_ren/CTCF/01.QC/01.analysis/CTCF-control/CTCF-control.CPM.bw --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/integrate/integrate.yaml -p --allowed-rules atac_bw --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.aq3pmbge/14.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.aq3pmbge/14.jobfailed; exit 1)

