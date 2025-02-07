#!/usr/bin/sh
#snakemake -ps /work/frasergen/3D/work/shaojie/script/HiC/integrate/integrate_snake.py --configfile multi_sample.yaml --cores 1 #--cluster 'sbatch -p xhacexclu03,xhacexclu12,xhacexclu16 -n 1 -N 1 --cpus-per-task=1' --jobs 1

snakemake -ps pgt_snake.py --configfile pgt.yaml --cores 24 --cluster "sbatch -p xhacexclu12,xhacexclu16,xhacexclu03 -n 1 -c 1" --jobs 24 --rerun-incomplete

###perl -anle 'BEGIN{open(IN,"<","target.genelist.bak");%a=map{($i)=split " ";$i=>1}<IN>;close IN} $F[3]=~/.*(?=.ITAG)/;print $F[3] if $a{$&}' /work/frasergen/backup/3d/project/Interactome/231492_fanqie_20231110/04.Hi-C/01.ref/genome_gene.bed > target.genelist

#snakemake -ps integrate_snake.py --configfile integrate.yaml -j 10 --cluster 'sbatch -p xhacexclu03,xhacexclu12,xhacexclu16 -n 1 -N 1 --cpus-per-task=1' --jobs 10


