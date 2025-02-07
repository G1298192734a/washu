#!/usr/bin/sh
#snakemake -ps hic_pre.snake.py --configfile pre.yaml -j 1

snakemake -ps integrate_snake.py --configfile integrate.yaml -j 1

