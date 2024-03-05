# 2bit和genomesize
faToTwoBit /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/00.ref/Slividans_1326.fasta genome.2bit
twoBitInfo genome.2bit genome.size

# gene annotation
ln -s /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/00.ref/Slividans_1326.gff genome
python3 format_gff.py genome
/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/bgzip genome.refbed
/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/tabix -p bed genome.refbed.gz

# gc bw
hgGcPercent -win=1000 -file=genome.GC -wigOut genome genome.2bit -noDots
wigToBigWig genome.GC genome.size genome.gc.bigwig

# js和json
python3 washu_utils.py --ref Slividans_1326 annot_tracks 
python3 washu_utils.py --ref Slividans_1326 chromsize genome.size
python3 washu_utils.py --ref Slividans_1326 genome_js whfs-xs-232291
python3 washu_utils.py --ref Slividans_1326 hub_json whfs-xs-232291

