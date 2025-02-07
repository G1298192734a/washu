#!/bin/sh
# properties = {"type": "single", "rule": "get_gtf", "local": false, "input": [], "output": ["/work/frasergen/backup/3d/project/Interactome/232188_Gossypium_hirsutum/02.RNA_ATAC/analysis20240618_washu/00.longest_gtf/InteracGhirsutum_mRNA_genename.gtf"], "wildcards": {}, "params": {"perl": "/public/frasergen/3D/work/wanghao/myconda/AGAT/bin/perl", "longest_gtf": "/public/frasergen/3D/work/wanghao/myconda/AGAT/bin/agat_sp_keep_longest_isoform.pl", "gffread": "/public/frasergen/PUB/software/gffread/gffread-0.12.6/gffread"}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 2, "cluster": {}}
 cd /work/frasergen/3D/work/shaojie/script/HiC/washU/upload && \
/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3 \
-m snakemake /work/frasergen/backup/3d/project/Interactome/232188_Gossypium_hirsutum/02.RNA_ATAC/analysis20240618_washu/00.longest_gtf/InteracGhirsutum_mRNA_genename.gtf --snakefile /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/upload_snake.py \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.em0kyv9h --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/upload_hic3.yaml -p --allowed-rules get_gtf --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /public/frasergen/PUB/software/python/Python-3.7.10/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.em0kyv9h/2.jobfinished || (touch /work/frasergen/3D/work/shaojie/script/HiC/washU/upload/.snakemake/tmp.em0kyv9h/2.jobfailed; exit 1)

