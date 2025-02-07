#!/bin/sh
# properties = {"type": "single", "rule": "cpt_motif", "local": false, "input": ["/work/frasergen/3D/work/shaojie/script/HiC/integrate/analysis_HiC_ATAC_RNA/01.integrate/CapHiC-test.minus.CapHiC-control/HiC_CUT_H3K27ac/02.compartment_AB_switch_peak/work.sh"], "output": [], "wildcards": {"diff": "CapHiC-test.minus.CapHiC-control", "atac": "HiC_CUT_H3K27ac"}, "params": {}, "log": ["/work/frasergen/3D/work/shaojie/script/HiC/integrate/analysis_HiC_ATAC_RNA/01.integrate/CapHiC-test.minus.CapHiC-control/HiC_CUT_H3K27ac/03.compartment_AB_peak_motif/work.sh"], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 15, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/integrate && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3.7 \
-m snakemake /work/frasergen/3D/work/shaojie/script/HiC/integrate/analysis_HiC_ATAC_RNA/01.integrate/CapHiC-test.minus.CapHiC-control/HiC_CUT_H3K27ac/03.compartment_AB_peak_motif/work.sh --snakefile /work/frasergen/3D/work/shaojie/script/HiC/integrate/integrate_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.m6q7oqp8 /work/frasergen/3D/work/shaojie/script/HiC/integrate/analysis_HiC_ATAC_RNA/01.integrate/CapHiC-test.minus.CapHiC-control/HiC_CUT_H3K27ac/02.compartment_AB_switch_peak/work.sh --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/integrate/integrate.yaml -p --allowed-rules cpt_motif --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.m6q7oqp8/15.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/integrate/.snakemake/tmp.m6q7oqp8/15.jobfailed; exit 1)

