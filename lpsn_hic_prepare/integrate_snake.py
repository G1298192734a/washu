import yaml,snakemake,re
import pandas as pd
from pathlib import Path

DIR,Gbed=(Path(config[i]).resolve() for i in ("Dir","Gene_bed"))
gene_bed=pd.read_table(Gbed,engine="c",header=None)
genes=gene_bed[3].to_list()
spl_rev={j:i for i,j in config["Spls"].items()}
minus2vs={".minus.".join(map(lambda x:spl_rev[x],i.split("-vs-"))):i for i in config["Spl_vs"].split(",")}
SPL,MINUS=map(list,(config["Spls"],minus2vs))

rule all:
  input:
    DIR.joinpath("01.genome_gene_express","work.sh"),
    expand(DIR.joinpath("02.diff_gene_express","{diff}.WashU.gene.LogFC"),diff=MINUS),
    expand("{depth}",depth=[DIR.joinpath("03.BAM_2_coverage",f"{i.stem}.SiteDepth.gz") for i in Path(config["RNA_bam"]).glob("*bam")]),
    expand(DIR.joinpath("03.BAM_2_coverage","{spl}.fill"),spl=SPL),

rule gbed:
  input:
    Gbed
  output:
    DIR.joinpath("genome_gene.bed")
  run:
    gene_bed.iloc[:,0:4].to_csv(output[0],sep="\t",header=None,index=None)

rule fpkm:
  input:
    config["Fpkm"],
    rules.gbed.output
  output:
    DIR.joinpath("01.genome_gene_express","work.sh")
  params:
    "/public/frasergen/3D/work/wanghao/pipeline/Hi-C_2.0_RNA_prepare/fpkm4washu_v1.0.py",
  run:
    for i in SPL:
      fpkm=pd.read_table(input[0],engine="c",usecols=lambda x:re.search(rf"gene.*id|{config['Spls'][i]}.*",x,re.IGNORECASE),index_col=[0])
      fpkm[fpkm.index.isin(genes)].mean(axis=1).to_csv(Path(*output).parent.joinpath(f"{i}.fpkm"),sep="\t",header=False)
      shell(f"""python3 {params[0]} --gene-bed {input[1]} --fpkm {Path(*output).parent.joinpath(f"{i}.fpkm")} --gene-fpkm {Path(*output).parent.joinpath(f"{i}.fpkm4washu")}""")
    script=rf"""#!/usr/bin/sh
for i in {" ".join(SPL)};do
  cut -f1,6-8 {input[0]}|awk 'NR == FNR{{genes[$4];next}}{{if($1 in genes){{total=0;for(i=2;i<=NF;i++){{total += $i}} print $1 "\t" total/(NF-1)}}}}' ../genome_gene.bed - > $i.fpkm && \
  python3 {params[0]} --gene-bed {input[1]} --fpkm $i.fpkm --gene-fpkm $i.fpkm4washu
done
"""
    open(*output,"w").write(script)

rule deg:
  input:
    lambda x:re.sub(r"[^./]+-vs-[^./]+",minus2vs.get(x.diff),config["DEG"]),
    rules.gbed.output
  output:
    DIR.joinpath("02.diff_gene_express","{diff}.gene.LogFC"),
    DIR.joinpath("02.diff_gene_express","{diff}.WashU.gene.LogFC"),
  params:
    "/public/frasergen/3D/work/wanghao/pipeline/Hi-C_2.0_RNA_prepare/fpkm4washu_v1.0.py",
    DIR.joinpath("02.diff_gene_express","work.sh")
  run:
    deg=pd.read_table(input[0],engine="c",usecols=lambda x:re.search(r"gene.*id|log.*F.*c",x,re.IGNORECASE),index_col=[0]).dropna()
    deg[deg.index.isin(genes)].to_csv(output[0],sep="\t",header=False)
    shell("python3 {params[0]} --gene-bed {input[1]} --fpkm {output[0]} --gene-fpkm {output[1]}")
    script=rf"""awk 'NR == FNR{{genes[$4];next}}{{if($1 in genes){{print $1 "\t" $5}}}}' ../genome_gene.bed {input[0]} > {wildcards.diff}.gene.LogFC && \
python3 {params[0]} --gene-bed ../genome_gene.bed --fpkm {wildcards.diff}.gene.LogFC --gene-fpkm {wildcards.diff}.WashU.gene.LogFC
#sort -k1,1 -k2,2n {wildcards.diff}.WashU.gene.LogFC > {wildcards.diff}.LogFC.sorted && /public/frasergen/PUB/software/htslib/htslib-1.12/bin/bgzip {wildcards.diff}.LogFC.sorted && /public/frasergen/PUB/software/htslib/htslib-1.12/bin/tabix -p bed {wildcards.diff}.LogFC.sorted.gz
"""
    if not Path(params[1]).is_file(): open(params[1],"w").write(script)

rule rna_bam:
  input:
    lambda x:Path(config["RNA_bam"]).joinpath(f"{Path(x.depth).stem[:-10]}.bam")
  output:
    "{depth}"
  params:
    "/work/frasergen/3D/work/shaojie/PanDepth/bin/pandepth"
  run:
    shell("{params} -i {input} -a -o %s" % DIR.joinpath("03.BAM_2_coverage",Path(*input).stem))

rule rna_depth:
  input:
    lambda x:[*DIR.glob(f"03.BAM_2_coverage/{config['Spls'][x.spl]}*Depth.gz")]
  output:
    DIR.joinpath("03.BAM_2_coverage","{spl}.fill"),
  params:
    DIR.joinpath("03.BAM_2_coverage","work.sh")
  run:
    depth=pd.concat([pd.read_table(i,engine="c",na_filter=False,low_memory=False,usecols=[2],header=None) for i in input],axis=1)
    depth=pd.DataFrame({"chr":["chr1"]*depth.shape[0],"site":[*depth.index[1:],depth.shape[0]],"depth":depth.mean(axis=1).values})
    depth.to_csv(*output,sep="\t",index=None,header=None)
    script=rf"""#!/usr/bin/sh
for i in {" ".join(config["Spls"][wildcards.spl])};do
  zcat $i-1.SiteDepth.gz $i-2.SiteDepth.gz|perl -anle 'next if !/chr/;$a[$F[1]+1]+=$F[2];END{{print join "\n",map{{join "\t","chr1",$_,($a[$_-1] ? $a[$_-1] : 0.0001)/2}} 1..$#a}}' > ${{i/RT/HT}}.fill
done
"""
    if not Path(*params).is_file(): open(*params,"w").write(script)


