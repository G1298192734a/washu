## cpg island(column name:https://genome.ucsc.edu/cgi-bin/hgTables?hgta_doSchemaDb=hg19&hgta_doSchemaTable=cpgIslandExt)
#faToTwoBit hg38.fa hg38_fa.2bit
#twoBitToFa hg38_fa.2bit stdout | maskOutFa stdin hard stdout \
#    | cpg_lh /dev/stdin 2> cpg_lh.err \
#    |  awk '{$2 = $2 - 1; width = $3 - $2;  printf("%s\t%d\t%s\t%s %s\t%s\t%s\t%0.0f\t%0.1f\t%s\t%s\n", $1, $2, $3, $5, $6, width, $6, width*$7*0.01, 100.0*2*$6/width, $7, $9);}' \
#    | sort -k1,1 -k2,2n > cpgIsland.bed

## cpg shore
#samtools faidx hg38.fa
#cut -f-2 hg38.fa.fai > hg38.fa.g
#bedtools slop -i cpgIsland.bed -g hg38.fa.g -b 2000 | bedtools merge \
#    | bedtools subtract -a - -b cpgIsland.bed > cpgShores.bed
## cpg shelve
#bedtools slop -i cpgIsland.bed -g hg38.fa.g -b 4000 | bedtools merge \
#    | bedtools subtract -a - -b cpgIsland.bed | bedtools subtract -a - -b cpgShores.bed > cpgShelves.bed

# repeat
#module load perl
#/public/frasergen/PUB/software/repeatmasker/RepeatMasker/RepeatMasker
#/work/frasergen/DNA/project/project_2024Q4/20240909_241683_220644_220429_seven_plants/06.plantLE02/03.annotation/01.repeat/denovo/repeatmasker/assembly_final.fasta.RepeatMasker.1107.shell
#/public/frasergen/PUB/software/repeatmasker/RepeatMasker-4.1.5/RepeatMasker
#目录下的famdb.py查找近源物种
/usr/bin/perl /public/frasergen/PUB/software/repeatmasker/RepeatMasker-4.1.2/RepeatMasker/RepeatMasker -nolow -no_is -norna -parallel 2 --species human chrY.fa


