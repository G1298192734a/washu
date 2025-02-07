#!/usr/bin/sh
#samtools sort -n -@ 16 Rdu-f0p-2.deduplicated.bam -o Rdu-f0p-2.deduplicated.sorted.bam
#samtools sort -n -@ 16 Tr01-14_B.deduplicated.bam -o Tr01-14_B.deduplicated.sorted.bam
#module load sambamba
#sambamba sort -m 8G -o Tr01-14_B.deduplicated.sorted1.bam Tr01-14_B.deduplicated.bam

#plotCoverage -b Tr01-14_B.deduplicated.sorted1.bam -p 10 -ct 0 -ct 1 -ct 2 -ct 5 -ct 10 -ct 15 -ct 20 -ct 30 -ct 50 --minMappingQuality 10 --smartLabels --samFlagExclude 256 --outCoverageMetrics genomeCoverage.coverageMetrics.txt -o genomeCoverage.txt > genomeCoverage.png
#plotCoverage -b Tr01-14_B.deduplicated.sorted1.bam -p 10 -o genomeCoverage.pdf
#pandepth -i Tr01-14_B.deduplicated.sorted1.bam -a -o Tr01-14_B

#MethylDackel mbias -@ 10 /work/frasergen/backup/3d/project/WGBS/241891_zhenjun/analysis_20241029/01.ref/ref_dir/genome.fasta Tr01-14_B.deduplicated.sorted1.bam Tr01-14_B

#MethylDackel extract -@ 10 /work/frasergen/backup/3d/project/WGBS/241891_zhenjun/analysis_20241029/01.ref/ref_dir/genome.fasta Tr01-14_B.deduplicated.sorted1.bam --opref Tr01-14_B

#/work/frasergen/backup/3d/project/WGBS/241891_zhenjun/analysis_20241029/01.ref/ref_dir/genome.gtf
#python3 gtf_extra_function_region.py /work/frasergen/backup/3d/project/WGBS/241891_zhenjun/analysis_20241029/01.ref/ref_dir/genome.gtf /work/frasergen/backup/3d/project/WGBS/241891_zhenjun/analysis_20241029/01.ref/ref_dir/genome.chromsizes
#computeMatrix reference-point --referencePoint center\
#  -R function_region.bed\
#  -S Tr01-14_B/Tr01-14_B.CG.bw Tr01-14_B/Tr01-14_B.CHG.bw Tr01-14_B/Tr01-14_B.CHH.bw\
#  -a 100 -b 100 -p 20 --skipZeros --missingDataAsZero \
#  -o Tr01-14_B_matrix.gz
computeMatrix reference-point --referencePoint center\
  -R genebody_region.bed\
  -S Tr01-14_B/Tr01-14_B.CG.bw Tr01-14_B/Tr01-14_B.CHG.bw Tr01-14_B/Tr01-14_B.CHH.bw\
  -a 100 -b 100 -p 20 --skipZeros --missingDataAsZero \
  -o Tr01-14_B_genebody_matrix.gz



