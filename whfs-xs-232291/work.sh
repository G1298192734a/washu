# 2bit和genomesize
#faToTwoBit /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/00.ref/Slividans_1326.fasta genome.2bit
#twoBitInfo genome.2bit genome.size

# gene annotation
#ln -s /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/00.ref/Slividans_1326.gff genome
#python3 format_gff.py genome
#/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/bgzip genome.refbed
#/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/tabix -p bed genome.refbed.gz

# gc bw
#hgGcPercent -win=1000 -file=genome.GC -wigOut genome genome.2bit -noDots
#wigToBigWig genome.GC genome.size genome.gc.bigwig

# heatmap
#/usr/bin/perl /public/frasergen/3D/work/shaojie/pipeline/ThreeDgenome/06.WashU_data_prepare/matrix2pairwise_filter.pl -m /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/05.HiC/04.map/Slividans_1326/bychr/reslu_5k/heatmap_bychr/corrected_chr1.txt -sta 0 -end 8496762 -chr chr1 -r 5000 -id Slividans_1326
#/public/frasergen/PUB/software/python/Python-2.7.18/bin/python /public/frasergen/3D/work/shaojie/pipeline/ThreeDgenome/06.WashU_data_prepare/pairwise2tabix.py Slividans_1326.chr1.0.8496762.matrix.pairwise Slividans_1326.5k.Heatmap
#sort -T ./tmp -k1,1 -k2,2n Slividans_1326.5k.Heatmap > Slividans_1326.5k.Heatmap.sorted
#/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/bgzip Slividans_1326.5k.Heatmap.sorted
#/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/tabix -p bed Slividans_1326.5k.Heatmap.sorted.gz

# compartment
#sort -k1,1 -k2,2n /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/05.HiC/15.compartment/Slividans_1326_100k/sample_100000.compartment_cscore.bedgraph.new > Slividans_1326.100k.Compartment.sorted
#/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/bgzip Slividans_1326.100k.Compartment.sorted
#/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/tabix -p bed Slividans_1326.100k.Compartment.sorted.gz

# tad
#awk '{print $1"\t"$2"\t"$3}' /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/05.HiC/06.tad/02.insulation/Slividans_1326_1k/TAD_genome.bed > Slividans_1326.1k.TAD1
#/public/frasergen/PUB/software/python/Python-2.7.18/bin/python /public/frasergen/3D/work/shaojie/pipeline/ThreeDgenome/06.WashU_data_prepare/bed2tabix.py Slividans_1326.1k.TAD1 Slividans_1326.1k.TAD
#sort -k1,1 -k2,2n Slividans_1326.1k.TAD > Slividans_1326.1k.TAD.sorted
#/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/bgzip Slividans_1326.1k.TAD.sorted
#/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/tabix -p bed Slividans_1326.1k.TAD.sorted.gz

# interaction
/usr/bin/perl /public/frasergen/3D/pipeline/Interactome/bacteria/hicare_v1/src/micro/script/interaction2pairwise.pl -i1 /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/05.HiC/05.loops/Slividans_1326_1k/02.annot/cis_interaction.sort -i2 /work/frasergen/backup/3d/project/Interactome/232291_lianmeijun_20240118/00.ref/Slividans_1326.fasta.fai -o1 Slividans_1326_1k_cis.pairwise -r 1000
/public/frasergen/PUB/software/python/Python-2.7.9/bin/python2.7 /public/frasergen/3D/pipeline/Interactome/bacteria/hicare_v1/src/micro/script/pairwise2tabix.py Slividans_1326_1k_cis.pairwise Slividans_1326_1k_cis.Interaction
sort -k1,1 -k2,2n Slividans_1326.1k.Interaction > Slividans_1326.1k.Interaction.sorted
/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/bgzip Slividans_1326.1k.Interaction.sorted
/public/frasergen/PUB/software/samtools/samtools-1.13/htslib-1.13/tabix -p bed Slividans_1326.1k.Interaction.sorted.gz

# js和json
#python3 washu_utils.py --ref Slividans_1326 annot_tracks 
#python3 washu_utils.py --ref Slividans_1326 chromsize genome.size
#python3 washu_utils.py --ref Slividans_1326 genome_js whfs-xs-232291
#python3 washu_utils.py --ref Slividans_1326 hub_json whfs-xs-232291

