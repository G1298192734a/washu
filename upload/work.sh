#!/usr/bin/sh
module load python/3.7.10 snakemake

#snakemake -ps upload_hic2.0_snake.py --configfile upload_hic2.yaml -j 4 --cluster 'sbatch -p xhacexclu12,xhacexclu16,xhacexclu03 -n 1 -N 1 --cpus-per-task=1' --jobs 4 --rerun-incomplete
#snakemake -ps upload_hic3.0_snake.py --configfile upload_hic3.yaml -j 4 --cluster 'sbatch -p xhacexclu12,xhacexclu16,xhacexclu03 -n 1 -N 1 --cpus-per-task=1' --jobs 4 --rerun-incomplete
snakemake -ps upload_hic3.0_snake.py --configfile upload_hic3.yaml -j 4

#snakemake -ps test.py --cores 1 -n

