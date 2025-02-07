import yaml,snakemake,re,shutil
import pandas as pd
import pyBigWig as pb
from pathlib import Path
from itertools import chain

DIR,Gbed=(Path(config[i]).resolve() for i in ("Dir","Gene_bed"))
gene_bed=pd.read_table(Gbed,engine="c",header=None)
genes=gene_bed[3].to_list()
PRE,INTEG=map(DIR.joinpath,("00.pre","01.integrate"))

diff_config = yaml.load(open(config["YAML"], "r"), Loader=yaml.FullLoader)
config = {**config, **diff_config}
HIC,DIFF=Gbed.parent.parent,Path(config["YAML"]).parent.parent
mapfile=dict((i.split() for i in open(HIC.parent.joinpath("00.ref","mapfile.txt"))))
minus2vs=dict(zip(*map(lambda x:config[x].split(","),("sample_minus","Spl_vs"))))
SPL,MINUS=map(list,(config["Spls"],minus2vs))
ATAC,atachash=(list(config["ATAC"]),{
  i:dict(zip(config["sample_minus"].split(","),j[1].split(",")))
  for i,j in config.get("ATAC").items()
}) if config.get("ATAC") else (list(),dict())

cpt_res,tad_res,loop_res=map(config.get,("compartment_res","tad_res","loop_res"))

rule all:
  input:
    expand(PRE.joinpath("01.genome_gene_express","{spl}.fpkm4washu"),spl=SPL),
    expand(PRE.joinpath("02.diff_gene_express","{diff}.WashU.gene.LogFC"),diff=MINUS),
    PRE.joinpath("03.BAM_2_coverage","work.sh"),
    expand(INTEG.joinpath("{diff}","HiC_RNA","01.compartment_RNA","Done"),diff=MINUS),
    expand(INTEG.joinpath("{diff}","HiC_RNA","02.TAD_RNA","Done"),diff=MINUS),
    expand(INTEG.joinpath("{diff}","{atac}","01.compartment_AB_peak","Done"),diff=MINUS,atac=ATAC),
    expand(INTEG.joinpath("{diff}","{atac}","02.compartment_AB_switch_peak","Done"),diff=MINUS,atac=ATAC),
    expand(INTEG.joinpath("{diff}","{atac}","03.compartment_AB_peak_motif","Done"),diff=MINUS,atac=ATAC),
    expand(INTEG.joinpath("{diff}","{atac}","04.tad_peak","Done"),diff=MINUS,atac=ATAC),
    expand(INTEG.joinpath("{diff}","{atac}","05.loop_peak","Done"),diff=MINUS,atac=ATAC),
    DIR.joinpath(re.sub("^analysis",f"whfs-xs-{HIC.parent.name[:6]}",config["Dir"]))

rule gbed:
  input:
    Gbed
  output:
    PRE.joinpath("genome_gene.bed")
  run:
    gene_bed.iloc[:,0:4].to_csv(output[0],sep="\t",header=None,index=None)

rule fpkm:
  input:
    config["Fpkm"],
    PRE.joinpath("genome_gene.bed")
  output:
    PRE.joinpath("01.genome_gene_express","{spl}.fpkm"),
    PRE.joinpath("01.genome_gene_express","{spl}.fpkm4washu"),
    PRE.joinpath("01.genome_gene_express","{spl}_biorep.txt")
  params:
    "/public/frasergen/3D/work/wanghao/pipeline/Hi-C_2.0_RNA_prepare/fpkm4washu_v1.0.py",
    PRE.joinpath("01.genome_gene_express","work.sh")
  run:
    fpkm=pd.read_table(input[0],engine="c",usecols=lambda x:re.search(r"gene.*id|{0}".format(config["Spls"][wildcards.spl]),x,re.IGNORECASE),index_col=[0])
    with open(output[2],"w") as REP:
      print(*(f"{wildcards.spl}\t{i}" for i in fpkm.columns),sep="\n",file=REP)
    fpkm[fpkm.index.isin(genes)].mean(axis=1).to_csv(output[0],sep="\t",header=False)
    script=rf"""cut -f1,6-8 {input[0]}|awk 'NR == FNR{{genes[$4];next}}{{if($1 in genes){{total=0;for(i=2;i<=NF;i++){{total += $i}} print $1 "\t" total/(NF-1)}}}}' ../genome_gene.bed - > {wildcards.spl}.fpkm && \
python3 {params[0]} --gene-bed ../genome_gene.bed --fpkm {wildcards.spl}.fpkm --gene-fpkm {wildcards.spl}.fpkm4washu
"""
    if not Path(params[1]).is_file(): open(params[1],"w").write(script)
    shell("python3 {params[0]} --gene-bed {input[1]} --fpkm {output[0]} --gene-fpkm {output[1]}")

def get_deg(diff,temple):
  return re.sub(r"[^./]+-vs-[^./]+",minus2vs.get(diff),temple)
rule deg:
  input:
    rules.gbed.output
  output:
    PRE.joinpath("02.diff_gene_express","{diff}.gene.LogFC"),
    PRE.joinpath("02.diff_gene_express","{diff}.WashU.gene.LogFC"),
    PRE.joinpath("02.diff_gene_express","{diff}.DEseq2.before-filter.xls")
  params:
    "/public/frasergen/3D/work/wanghao/pipeline/Hi-C_2.0_RNA_prepare/fpkm4washu_v1.0.py",
    PRE.joinpath("02.diff_gene_express","work.sh")
  run:
    Deg,prefilter=get_deg(wildcards.diff,config["DEG"]),get_deg(wildcards.diff,config["Filter_before"])
    deg=pd.read_table(Deg,engine="c",usecols=lambda x:re.search(r"gene.*id|log.*F.*c",x,re.IGNORECASE),index_col=[0]).dropna()
    deg[deg.index.isin(genes)].to_csv(output[0],sep="\t",header=False)
    script=rf"""awk -F"\t" 'NR == FNR{{genes[$4];next}}{{if($1 in genes){{print $1 "\t" $5}}}}' ../genome_gene.bed {Deg} > {wildcards.diff}.gene.LogFC && \
python3 {params[0]} --gene-bed ../genome_gene.bed --fpkm {wildcards.diff}.gene.LogFC --gene-fpkm {wildcards.diff}.WashU.gene.LogFC
fc=1
fdr=0.05
awk -v fc=$fc -v fdr=$fdr -v FS="\t" -v OFS="\t" 'NR==1{{for(i=1;i<=NF;i++){{if($i=="log2FoldChange" || $i=="logFC"){{col1=i}} else if($i=="pvalue" || $i=="Pvalue"){{col2=i}}}};print $0;next}}$col2!="NA"{{if(($col1>-fc && $col1<fc) || $col2>fdr) $col1="0.0";print$0}}' {prefilter}|sed 's/log2FoldChange/logFC/' > {wildcards.diff}.DEseq2.before-filter.xls
#sort -k1,1 -k2,2n {wildcards.diff}.WashU.gene.LogFC > {wildcards.diff}.LogFC.sorted && /public/frasergen/PUB/software/htslib/htslib-1.12/bin/bgzip {wildcards.diff}.LogFC.sorted && /public/frasergen/PUB/software/htslib/htslib-1.12/bin/tabix -p bed {wildcards.diff}.LogFC.sorted.gz
"""
    if not Path(params[1]).is_file(): open(params[1],"w").write(script)
    shell(r"""python3 {params[0]} --gene-bed {input[0]} --fpkm {output[0]} --gene-fpkm {output[1]} && \
awk -v fc=1 -v fdr=0.05 -v FS="\t" -v OFS="\t" 'NR==1{{for(i=1;i<=NF;i++){{if($i=="log2FoldChange" || $i=="logFC"){{col1=i}} else if($i=="pvalue" || $i=="Pvalue"){{col2=i}}}};print $0;next}}$col2!="NA"{{if(($col1>-fc && $col1<fc) || $col2>fdr) $col1="0.0";print$0}}' {prefilter}|sed 's/log2FoldChange/logFC/;s/gene_id/GeneID/' > {output[2]}
""")

rule rna_bw:
  output:
    PRE.joinpath("03.BAM_2_coverage","work.sh")
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
#python3 {params[1]} --insh script.sh nrows 4
"""
    open(output[0],"w").write(script)

def get_compt_gc(spl):
  return next(HIC.glob(f"*compart*/*PCA*/{spl}*{cpt_res}/*gc*bed"))
def get_cpt_switch(diff):
  return next(DIFF.glob(f"*compart*/{diff}/*switch.xls"))
rule compt_rna:
  input:
    PRE.joinpath("02.diff_gene_express","{diff}.DEseq2.before-filter.xls")
  output:
    PRE.joinpath("02.diff_gene_express","{diff}_biorep.txt"),
    INTEG.joinpath("{diff}","HiC_RNA","01.compartment_RNA","Done")
  log:
    INTEG.joinpath("{diff}","HiC_RNA","01.compartment_RNA","work.sh")
  params:
    "/public/frasergen/3D/pipeline/Interactome/multi_omics/HiC3.0/hic_rnaseq/compartment.r"
  run:
    treat,cont=wildcards.diff.split(".minus.")
    script=rf"""#!/usr/bin/sh
cat {PRE.joinpath("01.genome_gene_express",treat)}_biorep.txt {PRE.joinpath("01.genome_gene_express",cont)}_biorep.txt > {output[0]}
module load R/3.6.3 bedtools
Rscript {params[0]} --outpfix {cont}-vs-{treat} --splname1 {cont} --splname2 {treat} \
  --rna_biorep {output[0]} \
  --genebed {Gbed} \
  --compartment1 {get_compt_gc(cont)} \
  --compartment2 {get_compt_gc(treat)} \
  --switch {get_cpt_switch(wildcards.diff)} \
  --before_filter {input[0]} && touch Done
"""
    open(log[0],"w").write(script)
    shell(f"cd {Path(log[0]).parent} && bash work.sh")

def get_tad_gc(spl):
  return next(HIC.glob(f"*tad*/*cool*/{spl}*{tad_res}/*gc*density"))
rule tad_rna:
  input:
    PRE.joinpath("02.diff_gene_express","{diff}.DEseq2.before-filter.xls"),
    PRE.joinpath("02.diff_gene_express","{diff}_biorep.txt")
  output:
    INTEG.joinpath("{diff}","HiC_RNA","02.TAD_RNA","Done")
  log:
    INTEG.joinpath("{diff}","HiC_RNA","02.TAD_RNA","work.sh")
  params:
    "/public/frasergen/3D/pipeline/Interactome/multi_omics/HiC3.0/hic_rnaseq/tad.r"
  run:
    treat,cont=wildcards.diff.split(".minus.")
    script=rf"""#!/usr/bin/sh
module load R/3.6.3 bedtools
Rscript {params[0]} --outpfix {cont}-vs-{treat} --splname1 {cont} --splname2 {treat} \
  --rna_biorep {input[1]} \
  --genebed {Gbed} \
  --tad1 {get_tad_gc(cont)} \
  --tad2 {get_tad_gc(treat)} \
  --before_filter {input[0]} && touch Done
"""
    open(log[0],"w").write(script)
    shell(f"cd {Path(log[0]).parent} && bash work.sh")

def get_merg_bed(wildcards):
  return Path(config["ATAC"][wildcards.atac][0]).glob("00*/*peak/*bed")
rule atac_merg_bed:
  input:
    get_merg_bed
  output:
    directory(PRE.joinpath("{atac}_pre","00.merg_bed"))
  run:
    merg_bed=Path(output[0])
    merg_bed.mkdir(parents=True, exist_ok=True)
    for i in input:
      mergbed=pd.read_table(i,engine="c",header=None)
      mergbed[0]=mergbed[0].replace(mapfile)
      mergbed.to_csv(merg_bed.joinpath(Path(i).name),sep="\t",index=None,header=None)
  
def get_diff_peak(wildcards):
  return Path(config["ATAC"][wildcards.atac][0]).glob("*D*/*/*D*tsv")
rule atac_diff_tsv:
  input:
    get_diff_peak
  output:
    directory(PRE.joinpath("{atac}_pre","01.diff_tsv"))
  run:
    diff_tsv=Path(output[0])
    diff_tsv.mkdir(parents=True, exist_ok=True)
    for i in input:
      i=Path(i)
      if re.search("Gain|Loss",i.name):continue
      difftsv=pd.read_table(i,engine="c")
      difftsv["seqnames"]=difftsv["seqnames"].replace(mapfile)
      difftsv.to_csv(diff_tsv.joinpath(i.name),sep="\t",index=None)
  
def get_cpm_bw(wildcards):
  return Path(config["ATAC"][wildcards.atac][0]).parent.glob("01*/01*/*/*CPM.bw")
rule atac_bw:
  input:
    get_cpm_bw
  output:
    directory(PRE.joinpath("{atac}_pre","02.cpm_bw"))
  run:
    cpm_bw=Path(output[0])
    cpm_bw.mkdir(parents=True, exist_ok=True)
    for i in input:
      bw=pb.open(i)
      hdr=[(mapfile[chrom], length) for chrom, length in bw.chroms().items() if chrom in mapfile]
      bwOutput = pb.open(str(cpm_bw.joinpath(Path(i).name)), "w")
      bwOutput.addHeader(hdr)
      for chrom, length in bw.chroms().items():
        ints = bw.intervals(chrom, 0, length)
        if not (chrom in mapfile and len(ints)):continue
        bwOutput.addEntries([mapfile[chrom]] * len(ints), [x[0] for x in ints], ends=[x[1] for x in ints], values=[x[2] for x in ints])
      bw.close()
      bwOutput.close()

def get_input_bed(wildcards):
  return [PRE.joinpath(f"{wildcards.atac}_pre","00.merg_bed",f"{i}.merge_peak.bed") for i in atachash[wildcards.atac][wildcards.diff].split("-vs-")]
rule cpt_atac:
  input:
    rules.atac_merg_bed.output
  output:
    INTEG.joinpath("{diff}","{atac}","01.compartment_AB_peak","Done")
  log:
    INTEG.joinpath("{diff}","{atac}","01.compartment_AB_peak","work.sh")
  params:
    "/public/frasergen/3D/pipeline/Interactome/multi_omics/HiC3.0/hic_peak/cptmt_peak.r"
  run:
    treat,cont=wildcards.diff.split(".minus.")
    peaks=get_input_bed(wildcards)
    script=rf"""#!/usr/bin/sh
module load R/3.6.3 bedtools
Rscript {params[0]} --outdir ./ --samplename {treat} \
  --compartment {get_compt_gc(treat)} \
  --peak {peaks[1]} && \
Rscript {params[0]} --outdir ./ --samplename {cont} \
  --compartment {get_compt_gc(cont)} \
  --peak {peaks[0]} && touch Done
"""
    open(log[0],"w").write(script)
    shell(f"cd {Path(log[0]).parent} && bash work.sh")

def get_input_tsv(wildcards):
  return next(PRE.glob(f"{wildcards.atac}*/01*/{atachash[wildcards.atac][wildcards.diff]}*"))
rule cpt_switch_atac:
  input:
    rules.atac_diff_tsv.output
  output:
    INTEG.joinpath("{diff}","{atac}","02.compartment_AB_switch_peak","Done")
  log:
    INTEG.joinpath("{diff}","{atac}","02.compartment_AB_switch_peak","work.sh")
  params:
    "/work/frasergen/3D/pipeline/Interactome/multi_omics/HiC3.0/hic_peak/cptmt_diffpeak.r"
  run:
    treat,cont=wildcards.diff.split(".minus.")
    script=rf"""#!/usr/bin/sh
module load R/3.6.3 bedtools
Rscript {params[0]} --outpfix {atachash[wildcards.atac][wildcards.diff]} --splname1 {cont} --splname2 {treat} \
  --compartment1 {get_compt_gc(cont)} \
  --compartment2 {get_compt_gc(treat)} \
  --switch {get_cpt_switch(wildcards.diff)} \
  --diff_peak {get_input_tsv(wildcards)} && touch Done
"""
    open(log[0],"w").write(script)
    shell(f"cd {Path(log[0]).parent} && bash work.sh")

def get_merg_fa(atac):
  return next(Path(config["ATAC"][atac][0]).glob("00*/*peak.fa"))
rule cpt_motif:
  input:
    INTEG.joinpath("{diff}","{atac}","02.compartment_AB_switch_peak","Done")
  output:
    INTEG.joinpath("{diff}","{atac}","03.compartment_AB_peak_motif","Done")
  log:
    INTEG.joinpath("{diff}","{atac}","03.compartment_AB_peak_motif","work.sh")
  params:
    "/work/frasergen/3D/pipeline/Interactome/multi_omics/HiC3.0/hic_peak/generate_motif.py"
  run:
    treat,cont=wildcards.diff.split(".minus.")
    adiff=atachash[wildcards.atac][wildcards.diff]
    script=rf"""#!/usr/bin/sh
jaspar=$(perl -nle 'next if not /(?<=^MEME: ).*/;print $&' {Path(config["ATAC"][wildcards.atac][0]).parent.joinpath("project.yaml")})
awk 'NR==1{{print $0;next}}$4 == "A2B"|| $4 =="B2A"' ../02.compartment_AB_switch_peak/*switch_DiffPeak.xls > {wildcards.diff}_switch_DiffPeak.xls && \
/public/frasergen/PUB/software/python/Python-3.7.7/bin/python3 {params[0]} \
  --outdir ./ --diff {adiff} --switch_diffpeak {wildcards.diff}_switch_DiffPeak.xls \
  --inref {HIC.joinpath("01.ref","genome.fa")} \
  --control_fa {get_merg_fa(wildcards.atac)} \
  --jaspar $jaspar && touch Done
"""
    open(log[0],"w").write(script)
    shell(f"cd {Path(log[0]).parent} && bash work.sh")

rule tad_peak:
  input:
    rules.atac_bw.output
  output:
    INTEG.joinpath("{diff}","{atac}","04.tad_peak","Done")
  log:
    INTEG.joinpath("{diff}","{atac}","04.tad_peak","work.sh")
  params:
    "/work/frasergen/3D/pipeline/Interactome/multi_omics/HiC3.0/hic_peak/tad_diffpeak.py"
  run:
    treat,cont=wildcards.diff.split(".minus.")
    adiff=atachash[wildcards.atac][wildcards.diff]
    agroup=pd.read_csv(Path(config["ATAC"][wildcards.atac][0]).parent.joinpath("samplelist"),sep="\t",engine="c",header=None,usecols=[0,1])
    agroup[0]=agroup[0].apply(lambda x:str(next(PRE.glob(f"{wildcards.atac}*/02*/{x}*bw"))))
    agroup=agroup.groupby(1)
    bws=[*map(lambda x:agroup.get_group(x)[0].str.cat(sep=" "),adiff.split("-vs-"))]
    script=rf"""#!/usr/bin/sh
t="{bws[1]}"
t_n=$(echo $t|sed 's/ /\n/g'|xargs -n1 -I [] basename [] .bw)
c="{bws[0]}"
c_n=$(echo $c|sed 's/ /\n/g'|xargs -n1 -I [] basename [] .bw)
module load deeptools
/public/frasergen/PUB/software/python/Python-3.7.7/bin/python3 {params[0]} \
  --outpfix {adiff}_tad --splname1 {cont} --splname2 {treat} \
  --tad1 {get_tad_gc(cont)} \
  --tad2 {get_tad_gc(treat)} \
  --bw1 $t --bw1_names $t_n --bw2 $c --bw2_names $c_n && touch Done
"""
    open(log[0],"w").write(script)
    shell(f"cd {Path(log[0]).parent} && bash work.sh")

def get_loop(spl):
  return next(HIC.glob(f"07*/02*/{spl}/*combine.mustache")) if not config.get("loop") or config["loop"]=="mustache" else next(HIC.glob(f"07*/03*/{spl}*{loop_res}/04*/*cis*gz"))
rule loop_peak:
  input:
    rules.atac_diff_tsv.output
  output:
    INTEG.joinpath("{diff}","{atac}","05.loop_peak","Done")
  log:
    INTEG.joinpath("{diff}","{atac}","05.loop_peak","work.sh")
  params:
    "/work/frasergen/3D/pipeline/Interactome/multi_omics/HiC3.0/hic_peak/loop_diffpeak.py"
  run:
    treat,cont=wildcards.diff.split(".minus.")
    adiff=atachash[wildcards.atac][wildcards.diff]
    get_bed="\n".join((
      rf"""zcat {i} | awk 'BEGIN{{res={loop_res}}}{{a1_s=$2 - res/2 ;a1_e=$2 + res/2;a2_s=$4 - res/2 ;a2_e=$4 + res/2;print "chr" $1 "\t" a1_s  "\t" a1_e "\t" "chr" $3 "\t" a2_s  "\t" a2_e}}' > {j}_cis.bedpe""" if i.name.endswith("gz") else f"sed 1d {i}|cut -f1-6 > {j}_cis.bedpe"
      for i,j in zip(map(get_loop,(treat,cont)),(treat,cont))
))
    script=rf"""#!/usr/bin/sh
{get_bed}
/public/frasergen/PUB/software/python/Python-3.7.7/bin/python3 {params[0]} \
  --genebed {Gbed} \
  --splname1 {cont} --splname2 {treat} --outpfix {adiff}_loop_dar \
  --loop1 {cont}_cis.bedpe --loop2 {treat}_cis.bedpe \
  --diff_peak {get_input_tsv(wildcards)} && touch Done
"""
    open(log[0],"w").write(script)
    shell(f"cd {Path(log[0]).parent} && bash work.sh")

rule result:
  input:
    expand(INTEG.joinpath("{diff}","HiC_RNA","01.compartment_RNA","Done"),diff=MINUS),
    expand(INTEG.joinpath("{diff}","HiC_RNA","02.TAD_RNA","Done"),diff=MINUS),
    expand(INTEG.joinpath("{diff}","{atac}","01.compartment_AB_peak","Done"),diff=MINUS,atac=ATAC),
    expand(INTEG.joinpath("{diff}","{atac}","02.compartment_AB_switch_peak","Done"),diff=MINUS,atac=ATAC),
    expand(INTEG.joinpath("{diff}","{atac}","03.compartment_AB_peak_motif","Done"),diff=MINUS,atac=ATAC),
    expand(INTEG.joinpath("{diff}","{atac}","04.tad_peak","Done"),diff=MINUS,atac=ATAC),
    expand(INTEG.joinpath("{diff}","{atac}","05.loop_peak","Done"),diff=MINUS,atac=ATAC),
  output:
    directory(DIR.joinpath(re.sub("^analysis",f"whfs-xs-{HIC.parent.name[:6]}",config["Dir"])))
  params:
    "/work/frasergen/3D/work/shaojie/script/HiC/integrate/Hi-C_RNA_联合分析说明.docx",
    "/work/frasergen/3D/work/shaojie/script/HiC/integrate/Hi-C_ATAC_联合分析说明.docx"
  run:
    shutil.copytree(INTEG,output[0])
    for i in Path(output[0]).glob("*"):
      i.rename(i.parent.joinpath("-vs-".join(i.name.split(".minus.")[::-1])))
    for i in Path(output[0]).rglob("*"):
      if i.is_dir(): continue
      if re.search("-vs-|heatmap|^ame",i.name): continue
      i.unlink()
    shutil.copy(params[0],output[0])
    if config.get("ATAC"): shutil.copy(params[1],output[0])
    shell("find {output} -type f -print0|xargs -0 md5sum|sort -k2,2 > {output}/md5.txt")



