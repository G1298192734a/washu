import yaml,snakemake,re
from pathlib import Path
from itertools import chain
from configparser import ConfigParser

DIR,HIC,FAI=(Path(config[i]).resolve() for i in ("Dir","HiC","Fai"))
PRE,PGT=map(DIR.joinpath,("00.prepare_data","01.PGT"))

diff_config = yaml.load(open(config["YAML"], "r"), Loader=yaml.FullLoader)
config = {**config, **diff_config}
DIFF=Path(config["YAML"]).parent.parent
MINUS=config["sample_minus"].split(",")
GROUP=set(chain(*(i.split(".minus.") for i in MINUS)))
CHROM=[i.split()[0] for i in open(FAI)]
gen_res,cpt_res,tad_res,loop_res=map(config.get,("genome_res","compartment_res","tad_res","loop_res"))

def get_vs(minus):
  return "_vs_".join(minus.split(".minus.")[::-1])
def get_tad_cool(wildcards):
  return [*HIC.glob(f"06*/01*/{wildcards.spl}*{tad_res}*cool")][0]
def get_tad_bed(wildcards):
  return [*HIC.glob(f"06*/02*/{wildcards.spl}*{tad_res}/TAD_genome.bed")][0]
def get_tad_bw(wildcards):
  return [*HIC.glob(f"06*/02*/*/{wildcards.spl}*{tad_res}*bw")][0]
def get_cpt_bed(wildcards):
  return [*HIC.glob(f"05*/02*/*/{wildcards.spl}*{cpt_res}*compart*bed")][0]
def get_loop(wildcards):
  return [*HIC.glob(f"07*/02*/{wildcards.spl}/*combine.mustache")][0] if not config.get("loop") or config.get("loop")=="mustache" else [*HIC.glob(f"07*/03*/{wildcards.spl}*{loop_res}/04*/*cis*gz")][0]
def get_diff_mat(wildcards):
  return [*DIFF.glob(f"01*/02*/*/*/{'.minus.'.join(wildcards.diff.split('_vs_')[::-1])}*{wildcards.chr}.matrix.gz")][0]

rule all:
  input:
    expand(PRE.joinpath("01.cool_bychr","{spl}_"+f"{tad_res//1000}k.balance.cool"),spl=GROUP),
    expand(PRE.joinpath("02.minus_matrix_cool","{diff}","{chr}_zscore_minus.cool"),diff=[*map(get_vs,MINUS)],chr=CHROM),
    expand(PRE.joinpath("03.compartment","{spl}_"+f"{cpt_res//1000}k_cpt.bed"),spl=GROUP),
    expand(PRE.joinpath("04.tad","{diff}","corrected_{chr}.merged.TAD.bed"),diff=[*map(get_vs,MINUS)],chr=CHROM),
    expand(PRE.joinpath("05.tad_scores","{spl}."+f"{tad_res}.score.bw"),spl=GROUP),
    expand(PRE.joinpath("06.loop","{spl}_cis.bedpe"),spl=GROUP),
    PGT.joinpath("step01.target_region_plots.sh"),
    PRE.joinpath("07.bw","work.sh")

rule pgt:
  output:
    PGT.joinpath("step01.target_region_plots.sh"),
    PGT.joinpath("bigwigs.ini"),
  params:
    HIC.joinpath("01.ref","genome_gene.bed"),
    "/work/frasergen/3D/work/shaojie/script/HiC/integrate/target_genes_TAD_PGT_draw.sh",
    "/work/frasergen/3D/work/shaojie/script/HiC/integrate/bigwigs.ini"
  run:
    if config.get("BW"):
      for i,j in config["BW"].items():
        parser=ConfigParser()
        parser.read(params[2])
        parser.set("treat bw","file",j[0])
        parser.set("treat bw","title",f"$treat {i}")
        parser.set("cont bw","file",j[1])
        parser.set("cont bw","title",f"$cont {i}")
        with open(output[1],"a") as OUT:
          parser.write(OUT)
    else:
      shell("touch {output[1]}")
    pgts="\n".join((f"sh target_genes_TAD_PGT_draw.sh $(pwd)/draw_region.txt {re.sub('_vs_',' ',get_vs(i))}" for i in MINUS))
    script=rf"""#!/usr/bin/sh
### make draw region list
#chr10   27700000        28400000        NA      target  region
#awk '{{print $1 "\t" $2 "\t" $3 "\tNA\t" $4 "\t" $1 "_" $2 "_" $3 }}' gene.bed > draw_region.txt
#awk -v OFS="\t" 'NR==FNR{{genes[$1];next}}$4 in genes{{print $1,$2,$3,"NA",$4,$1 "_" $2 "_" $3,"region"}}' target.genelist {params[0]} > draw_region.txt
ln -s {FAI} genome.fa.fai;ln -s {re.sub("fasta.fai","gtf",str(FAI))} genome.gtf
cp {params[1]} . && sed -i '160 r bigwigs.ini' target_genes_TAD_PGT_draw.sh
### now draw plots!
{pgts}
#cp -r PGT_draw ../whfs-xs-{re.sub("^analysis",HIC.parent.name[:6],config["Dir"])}
#ls ../whfs*/*/*|grep -vP "png|pdf"|xargs rm -r
"""
    open(output[0],"w").write(script)

rule single:
  input:
    get_tad_cool
  output:
    PRE.joinpath("01.cool_bychr","{spl}_"+f"{tad_res//1000}k.balance.cool")
  params:
    PRE.joinpath("01.cool_bychr","work.sh")
  run:
    shell("ln -s {input} {output}")
    script=f"""#!/usr/bin/sh
spl="{" ".join(GROUP)}"
for i in $pl;do
  if [ ! -d $i ];then mkdir $i;fi
  hicTransform --matrix $i"_40k.balance.cool" --outFileName $i"_chr1_40k_oe.cool" --method obs_exp_non_zero --chromosomes chr1
  hicConvertFormat -m $i".hic" --inputFormat hic --outputFormat cool -o $i".cool" --resolutions 40000
  hicCorrectMatrix diagnostic_plot -m $i"_40000.cool" -o $i/$i"_matrix.pdf" &> $i/$i"_matrix_mad_threshold.out"
  madscore=$(grep "mad threshold" $i/$i"_matrix_mad_threshold.out" | sed 's/INFO:hicexplorer.hicCorrectMatrix:mad threshold //g');
  upper=$(echo -3*$madscore | bc);
  echo $madscore " " $upper > $i/$i"_matrix_mad_threshold.values"
  thresholds=$(cat $i/$i"_matrix_mad_threshold.values");
  hicCorrectMatrix correct --filterThreshold $thresholds -m $i"_40000.cool" -o $i/$i"_ICE_40000.cool" --correctionMethod ICE -n 500
  hicCorrectMatrix correct -m $i"_40000.cool" -o $i/$i"__40000.cool" --correctionMethod KR
done
hicCompareMatrices --matrices cont_chr1_40k_oe.cool treat_chr1_40k_oe.cool --outFileName cont_vs_treat_chr1_40k_oe.cool
"""
    if not Path(params[0]).is_file(): open(params[0],"w").write(script)

rule minus:
  input:
    get_diff_mat
  output:
    PRE.joinpath("02.minus_matrix_cool","{diff}","{chr}.matrix"),
    PRE.joinpath("02.minus_matrix_cool","{diff}","{chr}_zscore_minus.cool"),
  params:
    "/work/frasergen/3D/work/wanghao/my_scripts/python/hiclib_heatmap_dense_2_cool/dense_2_cool.py",
    PRE.joinpath("02.minus_matrix_cool","work.sh")
  run:
    chrsize=Path(output[0]).parent.parent.joinpath(f"{wildcards.chr}.chrom.sizes")
    shell(r"""awk -v chrom={wildcards.chr} '$1 == chrom{{print $1 "\t" $2}}' {FAI} > {chrsize} && \
zcat {input} | awk -F'\t' 'NR>1{{for(i=2;i<NF;i++){{printf("%.8f ",$i==""?0:$i)}};printf("%.8f\n",$NF==""?0:$NF)}}' > {output[0]} && \
/public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin/python {params[0]} {output[0]} {output[1]} {chrsize} {gen_res} float
""")
    script=rf"""#!/usr/bin/sh
res={gen_res}
length={len(CHROM)+1}
minus_id="{" ".join(map(get_vs,MINUS))}"
ref_fai={FAI}
zcat_gz(){{
  echo {re.sub(wildcards.chr,"chr$2",re.sub(wildcards.diff,"$1",input[0]))}
}}
for chrom in $(seq 1 $length); do
  awk -v chrom=chr$chrom '$1 == chrom{{print $1 "\t" $2}}' $ref_fai > chr$chrom".chrom.sizes"
done
for i in $minus_id; do 
  if [ ! -d $i ]; then mkdir $i; fi
  cd $i
  for chrom in $(seq 1 $length); do
    zcat $(zcat_gz $i $chrom)|awk -F'\t' 'NR>1{{for(i=2;i<NF;i++){{printf("%.8f ",$i==""?0:$i)}};printf("%.8f\n",$NF==""?0:$NF)}}' > chr$chrom".matrix"
    /public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin/python {params[0]} chr$chrom".matrix" chr$chrom"_zscore_minus.cool" ../chr$chrom".chrom.sizes" $res float
  done
  cd ..
done
"""
    if not Path(params[1]).is_file(): open(params[1],"w").write(script)

rule compt:
  input:
    get_cpt_bed
  output:
    PRE.joinpath("03.compartment","{spl}_"+f"{cpt_res//1000}k_cpt.bed")
  params:
    "/work/frasergen/backup/3d/project/Interactome/220006_zebrafish/07.personal/batch06/02.analysis/analysis_HiC_aftersale_20240530/03_PGT/00.prepare_data/03.compartment/parse_state_file.awk",
    PRE.joinpath("03.compartment","work.sh")
  run:
    #shell("awk -f {params} {input}|sort -k1,1V -k2,2V > {output}")
    shell(r"""perl -alne '$F[3]= $F[4]eq"A" ? 1 : -1;print join "\t",@F' {input} > {output}""")
    script=rf"""#!/usr/bin/sh
spl="{" ".join(GROUP)}"
bed(){{
  echo {re.sub(wildcards.spl,"$1",input[0])}
}}
for i in $spl;do
  perl -alne '$F[3]= $F[4]eq"A" ? 1 : -1;print join "\t",@F' $(bed $i)|sort -k1,1V -k2,2V > $i"_100k_cpt.bed"
done
"""
    if not Path(params[1]).is_file(): open(params[1],"w").write(script)

rule tad:
  input:
    get_tad_bed
  output:
    PRE.joinpath("04.tad","{spl}","corrected_{chr}.TAD.bed")
  params:
    PRE.joinpath("04.tad","work.sh")
  run:
    shell(r"""awk -v chrom={wildcards.chr} '$1 == chrom {{print $1 "\t" $2 "\t" $3}}' {input} > {output}""")
    script=rf"""#!/usr/bin/sh
length={len(CHROM)+1}
spl="{" ".join(GROUP)}"
diff="{" ".join(map(get_vs,MINUS))}"
bed(){{
  echo {re.sub(wildcards.spl,"$1",input[0])}
}}
for i in $spl;do
  if [ ! -d $i ]; then mkdir $i; fi
  cd $i
  for j in $(seq 1 $length);do
    awk -v chrom=chr$j '$1 == chrom {{print $1 "\t" $2 "\t" $3}}' $(bed $i) > corrected_chr$j".TAD.bed"
  done
  cd ..
done
for i in $diff;do
  if [ ! -d $i ]; then mkdir merged_tad_$i; fi
  spl1=${{i%%_vs_*}};spl2=${{i##*_vs_}}
  for j in $(seq 1 $length);do
    awk '!a[$0]++' $spl1/corrected_chr$j".TAD.bed" $spl2/corrected_chr$j".TAD.bed" > $i/corrected_chr$j".merged.TAD.bed"
  done
done
"""
    if not Path(params[0]).is_file(): open(params[0],"w").write(script)

rule merge_tad:
  input:
    expand(PRE.joinpath("04.tad","{spl}","corrected_{{chr}}.TAD.bed"),spl=GROUP)
  output:
    PRE.joinpath("04.tad","{diff}","corrected_{chr}.merged.TAD.bed")
  run:
    beds=[PRE.joinpath("04.tad",i,f"corrected_{wildcards.chr}.TAD.bed") for i in wildcards.diff.split("_vs_")]
    shell("awk '!a[$0]++' {beds} > {output}")

rule tad_bw:
  input:
    get_tad_bw
  output:
    PRE.joinpath("05.tad_scores","{spl}."+f"{tad_res}.score.bw")
  run:
    shell("ln -s {input} {output}")

rule loop:
  input:
    get_loop
  output:
    PRE.joinpath("06.loop","{spl}_cis.bedpe")
  params:
    PRE.joinpath("06.loop","work.sh"),
    "/public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin/pigz",
    "/public/frasergen/PUB/software/bedtools/bedtools-2.30.0/bin/bedtools",
  run:
    if input[0].endswith("gz"):
      #shell(r"""zcat {input} | awk 'BEGIN{{res={loop_res}}}{{a1_s=$2 - res/2 ;a1_e=$2 + res/2;a2_s=$4 - res/2 ;a2_e=$4 + res/2;print "chr" $1 "\t" a1_s  "\t" a1_e "\t" "chr" $3 "\t" a2_s  "\t" a2_e "\t" log($5)/log(10)}}' > {output}""")
      shell(r"""{params[1]} -p 30 -cd {input} | awk -v OFS="\t" '{{print "chr"$1,$2-half_res,$2+half_res,"chr"$3,$4-half_res,$4+half_res,log($5)/log(10)}}' half_res=$[{loop_res}/2] | {params[2]} sort -faidx {FAI} -i stdin > {output}""")
    else:
      shell("sed 1d {input} > {output}")
    script=rf"""#!/usr/bin/sh
res={loop_res}
spl="{" ".join(GROUP)}"
input(){{
  echo {re.sub(wildcards.spl,"$1",input[0])}
}}
for i in $spl;do
  {params[1]} -p 30 -cd $(input $i) | awk -v OFS="\t" '{{print "chr"$1,$2-half_res,$2+half_res,"chr"$3,$4-half_res,$4+half_res,log($5)/log(10)}}' half_res=$[$res/2] | {params[2]} sort -faidx {FAI} -i stdin > $i"_cis.bedpe"
  #sed 1d $(input $i) > $i"_cis.bedpe"
done
"""
    if not Path(params[0]).is_file(): open(params[0],"w").write(script)

rule bw:
  log:
    PRE.joinpath("07.bw","work.sh")
  params:
    "/public/frasergen/PUB/software/python/Python-3.7.10/bin/bamCoverage",
    "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/get_subscript.py"
  run:
    script=rf"""#!/bin/bash
for i in */*bam;do
  j=$(basename $i .bam)
  cat <<END >>script.sh
module load samtools/1.13 python/3.7.10
samtools view -H $i|sed -e 's#Superscaffold#chr#g' |samtools reheader - $i > $j.rename.bam && samtools index $j.rename.bam && \
{params[0]} -p 10 -b $j.rename.bam -bs 50 --normalizeUsing RPKM -o $j.bw && rm $j.rename.bam $j.rename.bam.bai
END
done
#for i in *bw;do
#j=$(basename $i .bw)
#if [ ! -d rename_bw ];then mkdir rename_bw;fi
#cat <<END >>script.sh
#bigWigToBedGraph $i stdout|perl -anle '\$F[0]=~/\\d+$/;\$chr=int($&)+1;\$F[0]="chr\$chr";print join "\\t",@F'|sort -k1,1 -k2,2n > $j.bdg && bedGraphToBigWig $j.bdg {FAI} rename_bw/$j.bw && rm $j.bdg
#END
#done
#python3 {params[1]} --insh script.sh nrows 4
"""
    open(*output,"w").write(script)



