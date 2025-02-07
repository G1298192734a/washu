#!/usr/bin/sh
module load python/3.7.10 snakemake

#snakemake -ps /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_snake.py --configfile trid1.yaml --cores 24 --cluster 'sbatch -p xhacexclu12,xhacexclu16,xhacexclu03 -n 1 -N 1 --cpus-per-task=1' --jobs 24 --rerun-incomplete

snakemake -ps /work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_snake.py --configfile trid.yaml --cores 24 --cluster 'sbatch -p xhacexclu16,xhacexclu13,xhacexclu12,xhacexclu20 -n 1 -N 1 --cpus-per-task=1' --jobs 24 --rerun-incomplete

#python3 get_marked.py /work/frasergen/3D/work/shaojie/script/HiC/3d_model/lianmeijun/02.3d_model/Slividans_1326/MDS.Slividans_1326.fa.fai-5000.pdb /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/05.HiC/01.ref/Slividans_1326.fa.fai /work/frasergen/3D/work/shaojie/script/HiC/3d_model/lianmeijun/01.matrix_fill/bin.bed 5000

