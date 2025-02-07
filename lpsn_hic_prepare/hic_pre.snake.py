#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os,time,yaml
from pathlib import Path

DIR=Path(config.get("DIR")).resolve()
project=f"whfs-xs-{os.path.basename(DIR).split('_',1)[0]}"
dirnames=("00.ref","01.rawdata","02.qc","03.enrich","04.RNA","05.HiC",)
REF,RAW,QC,ENRICH,RNA,HiC=map(lambda x:os.path.join(DIR,x),dirnames)
GENOME,MINUS=map(config.get,("GENOME","MINUS"))
SAMPLE,DATANAME=map(lambda x:config.get(x).split(","),("SAMPLE","DATANAME"))
group=[j for i in MINUS.split(",") for j in i.split("_vs_")]
GROUP=sorted(set(group),key=group.index) if group else SAMPLE
step=len(SAMPLE)//len(GROUP)
group_order="\n".join((f"  - {i}" for i in GROUP))

RNA=Path(RNA)
dirs=("01.genome_gene_express","02.diff_gene_express","03.BAM_2_coverage")
fpkm,deg,bam=map(lambda x:RNA.joinpath(x),dirs)

soft_config = yaml.load(open(config["PRE_SOFT_YAML"], "r"), Loader=yaml.FullLoader)
config = {**config, **soft_config}

date=time.strftime("%Y%m%d")

rule all:
  input:
    os.path.join(REF,"work.sh"),
    os.path.join(RAW,"work.sh"),
    expand(os.path.join(QC,"{sample}_qc.sh"),sample=SAMPLE),
    os.path.join(RNA,"work.sh"),
    os.path.join(ENRICH,"work.sh"),
#    os.path.join(HiC,"work.sh"),
    os.path.join(HiC,"report.sh"),
    os.path.join(HiC,"conf.yaml"),
    os.path.join(HiC,"report.yaml"),
#    os.path.join(REF,f"{GENOME}.mrna.fasta"),
#    expand(os.path.join(QC,"{sample}_trimming.sh"),sample=SAMPLE),
#    os.path.join(ENRICH,"prepare","enrich_prepare.sh"),

rule ref:
  output:
    os.path.join(REF,"work.sh"),
  params:
    get_map=config.get("get_map"),
    hic_util=config.get("hic_util"),
  run:
    script=rf"""##!/usr/bin/bash

## 准备参考和注释，细菌都只有一条染色体，需要判断是线性还是环状
# fa= #raw fasta
#gff= #raw gff
#chr1=
#g={GENOME}
#u={params.hic_util}

## 生成 mapfile
echo -e $chr1"\tchr1" > mapfile.txt

## fa/gff 取染色体部分 重命名 gene.bed gtf genemark mrna
#export PATH=/public/frasergen/PUB/software/python/Python-3.7.7/bin:$PATH
#module load gffread/0.12.6
#python3 $u fa_clean --infa $fa --chrom-map mapfile.txt --outpfix $g
#python3 $u gff_clean --ingff $gff --chrom-map mapfile.txt --outpfix $g && gffread -T $g.gff -o $g.gtf
#python3 $u gff_gene --ingff $g.gff --outpfix $g
#perl -i -anle 'print join "\t",@F[0..3]' ${{g}}_gene.bed
#python3 $u gtf_gene2mark --ingtf $g.gtf --outpfix $g
#gffread $g.gtf -g $g.fasta -x $g.mrna.fasta
    """
    open(*output,"w").write(script)

rule ref_sh:
  input:
    os.path.join(REF,"work.sh"),
  output:
    os.path.join(REF,f"{GENOME}.mrna.fasta"),
  shell:
    f"cd {REF} && bash work.sh"

rule rawdata:
  output:
    os.path.join(RAW,"work.sh"),
  params:
    get_rawdata=config.get("get_rawdata"),
    get_sub=config.get("get_sub"),
    data_storage=config.get("data_storage"),
  run:
    main=f"{project}_HiC_submit_{date}"
    sub="{00.rawdata,01.report,02.result}"
    
    f=lambda x:os.path.join(params.data_storage,project,main,*x)
    data_dir,submit_dir,data=map(f,(
      ["00.rawdata"], [sub], ["00*","*gz"],
    ))
    script=f"""##!/usr/bin/sh

## 准备数据目录及路径
#x=数据源导出.xlsx #默认从系统中下载数据文档
#mkdir -p {submit_dir}
#python3 {params.get_rawdata} $x

## 拆分数据合并脚本
#python3 {params.get_sub} --insh cp.sh numsign

## 生成备份数据 md5
##cat *md5.txt > md5.txt && rm *_md5.txt
##mv md5.txt *gz {data_dir} && ln -s {data} .

## 验证原始数据 md5,需要在合并完成后再次投递
##md5sum -c md5sum.txt > check.txt

    """
    open(*output,"w").write(script)

rule rawdata_sh:
  input:
    os.path.join(RAW,"work.sh"),
  output:
    os.path.join(RAW,"{sample}_2.fq.gz"),
  shell:
    f"cd {RAW} && bash work.sh"

rule qc:
  output:
    os.path.join(QC,"{sample}_qc.sh"),
  params:
    get_trim=config.get("get_trim"),
    figure_util=config.get("figure_util"),
    read1=os.path.join(RAW,"{sample}_R1.fq.gz"),
    read2=os.path.join(RAW,"{sample}_R2.fq.gz"),
  run:
    s=wildcards.sample
    figs=("per_base_quality.png","per_base_sequence_content.png","per_sequence_gc_content.png")
    f1=lambda x:os.path.join(s,f"{s}_clean_1P_fastqc","Images",x)
    f2=lambda x:os.path.join(s,f"{s}_clean_2P_fastqc","Images",x)
    files=" ".join([*map(f1,figs),*map(f2,figs)])
    script=rf"""#!/usr/bin/sh
export PATH=/public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin:$PATH
python3 {params.get_trim} \
  -sample {s} \
  -fq1 {params.read1} \
  -fq2 {params.read2} \
  -platform BGI \
  -out . \
  -threads 3
<<END
export PATH=/public/frasergen/PUB/software/python/Python-3.7.7/bin:$PATH
python3 {params.figure_util} combine_figure \
  --out-image {s}/{s}_clean/{s}.sequencing_quality.pdf \
  --ncol 3 {files}
python3 {params.figure_util} combine_figure \
  --out-image {s}/{s}_clean/{s}.sequencing_quality.png \
  --ncol 3 {files}
END
    """
    open(*output,"w").write(script)

rule qc_sh:
  input:
    script=os.path.join(QC,"{sample}_qc.sh"),
    read2=os.path.join(RAW,"{sample}_2.fq.gz"),
  output:
    os.path.join(QC,"{sample}_trimming.sh"),
  shell:
    "cd {QC} && bash {input.script}"

rule RNA:
  log:
    os.path.join(RNA,"work.sh")
  params:
    fpkm_pre=config["fpkm_pre"],
    fpkm=config["fpkm"],
    logfc=config["logfc"],
    integrate="/work/frasergen/3D/work/shaojie/script/HiC/lpsn_hic_prepare/integrate_snake.py",
  run:
    script=f"""#!/usr/bin/sh
#module load python/3.7.10 snakemake
cat <<END >integrate.yaml
Dir: ./
Gene_bed: {HiC.joinpath("01.ref",GENOME)}.bed
Fpkm: whfs*/02*/gene.fpkm_matrix*
DEG: whfs*/06*/genes/*known.DEG*
Spl_vs: rna_cont-vs-treat 
Spls:
{{0}}
YAML: {DIFF.joinpath("05.multi_report","report.yaml")}
END
snakemake -ps {params.integrate} --configfile integrate.yaml -j 10 --cluster 'sbatch -p xhacexclu03,xhacexclu12,xhacexclu16 -n 1 -c 1' --jobs 10
    """.format("\n".join((f"  {i}: {i}" for i in GROUP)))
    open(log[0],"w").write(script)

rule enrich:
  output:
    os.path.join(ENRICH,"work.sh"),
  params:
    enrich_prepare=config.get("enrich_prepare"),
    gene_mark=os.path.join(REF,f"{GENOME}.gene2mark"),
    mrna=os.path.join(REF,f"{GENOME}.mrna.fasta"),
  run:
    script=rf"""#!/usr/bin/bash
export PATH=/public/frasergen/PUB/software/Anaconda/anaconda3-3d/envs/iced/bin:$PATH
python {params.enrich_prepare} \
  -gene2mark {params.gene_mark} \
  -mfa {params.mrna} \
  -species b
    """
    open(*output,"w").write(script)

rule enrich_sh:
  input:
    os.path.join(REF,f"{GENOME}.mrna.fasta"),
  output:
    os.path.join(ENRICH,"prepare","enrich_prepare.sh"),
  shell:
    "cd {ENRICH} && bash work.sh"

rule HiC:
  output:
    hic=os.path.join(HiC,"work.sh"),
    report=os.path.join(HiC,"report.sh"),
  params:
    prepare=config.get("prepare_SG"),
    backup_bacteria=config.get("backup_bacteria"),
    micro_report=config.get("micro_report"),
    rmarkdown=config.get("rmarkdown"),
    moniter=config.get("moniter"),
    data_storage=os.path.join(config.get("data_storage"),project),
  run:
    script=f"""#!/usr/bin/bash
export PATH=/public/frasergen/PUB/software/python/Python-2.7.9/bin:$PATH
python {params.prepare} conf.yaml
#nohup perl {params.moniter} -list shell.list -o out.log -n 3 -t 60 -q xhacexclu12,xhacexclu03,xhacexclu16 &
    """
    open(output.hic,"w").write(script)
    script=f"""#!/usr/bin/bash
p={project}
d={date}

# make report result directory
export PATH=/public/frasergen/PUB/software/python/Python-2.7.9/bin:$PATH
python {params.backup_bacteria} report.yaml

# make rmd
module load python/3.7.10
python3 {params.micro_report} report report.yaml --outdir $p_rmd_report_$d

# make html
export R_LIBS=/public/frasergen/3D/work/wanghao/R_packages:/work/frasergen/PUB/software/R/R-3.6.3/lib64/R/library
module load R/3.6.3
Rscript {params.rmarkdown} -r ${{p}}_rmd_report_$d/$p.rmd -o $p.html

# dilivery {params.data_storage}
mv ${{p}}_rmd_report_$d/$p.rmd . && zip -r ${{p}}_HiC_report_$d.zip ${{p}}_rmd_report_$d
mv result ${{p}}_HiC_result_$d && zip -r ${{p}}_HiC_result_$d.zip ${{p}}_HiC_result_$d
"""
    open(output.report,"w").write(script)

rule hic_yaml:
  output:
    os.path.join(HiC,"conf.yaml"),
  params:
    motif_db=config.get("motif_db"),
  run:
    from itertools import chain
    rawdata=list()
    f=lambda x:(os.path.join(f"- {QC}",f"{x}_clean_{i}P.fastq.gz") for i in (1,2))
    for i,j in zip(range(0,len(SAMPLE),step),GROUP):
      a=((f"{p}:",*f(q)) for p,q in zip(SAMPLE[i:i+step],DATANAME[i:i+step]))
      rawdata+=(f"{j}:",*(f"  {p}" for p in chain.from_iterable(a)))
    rawdata="\n".join(rawdata)
    yaml=rf"""# reference: chromosome id must be sorted
# not include Mitochondria sequence and Chloroplast sequence
# just like this: >chr1 >chr2 >chr3 ...
reference: {os.path.join(REF,GENOME)}.fasta
chr_num: 1
sample_order:
{group_order}

{rawdata}

# step 0:all, 1:ref, 2:reads, 3:motif, 4:ice,
#    5:loops, 6:TAD, 7:correlate, 8:FIRE, 9:compact 
#    10:compare 11:enrich 12:fire_vs_rna-depth 13.washu 14.3D
step: 0
species: {GENOME}
threads: 6
outdir: .
enzyme: MboI
QC_tools: trimmomatic

# annot_bed: Chr-Start-End-GeneID (4column, header = F, sep = "\t")
annot_bed: {os.path.join(REF,GENOME)}_gene.bed

# fpkm_table: GeneID-Fpkm (2column, header = F, sep = "\t")
fpkm_table: 

# for step3 - motif analysis
motif_db: {params.motif_db}

# for step4 - ICE
# resolution unit: bp
# break unit: kb (<int>, >=1)
# map_method: map method choice from "ice" or "juicer"
genome_resolution: "1000,2000,5000"
bychr_resolution: "1000,2000,5000"
inter_distance_break: "200,400,600,800,1000,1200"
map_method: ice

# for step5 - call loops
# linear: genome is "True"(linear) or "False"(circle)
# target_gene: chrom-start-end-gene (4column, header = T, sep = "\t")
linear: False
loop_resolution: "1000,2000,5000"
loop_distance_break: "200,400,600,800,1000,1200"
#target_gene: 

# for step6 - call TAD
# bin_num: <int> (only for CID)
TAD_resolution: "1000,2000,5000"
bin_num: 10

# for step8 - FIRE(support 5000)
bin_range: 10000

# for step9 - compact(only support 5000)
compact_resolution: "5000"

# for step10 - compare
compare: {MINUS}
compare_genome_reslu: 5000
compare_bychr_reslu: 5000
compare_tad_method: insulation
compare_tad_reslu: 5000
compare_fithic_reslu: 1000
compare_compact_reslu: 5000
windows: 5
coefficient: 0.6
foldchange_table: 

# for step11 - enrich
GO: {os.path.join(ENRICH,"prepare","NR","we.go")}
KEGG: {os.path.join(ENRICH,"prepare","KEGG","kegg.ko")}

# for step12 - fire_vs_rna-depth
depth: 
"""
    open(*output,"w").write(yaml)

rule report_yaml:
  output:
    os.path.join(HiC,"report.yaml"),
  params:
    title=config.get("title"),
    species_info=config.get("species_info"),
    sample_info=config.get("sample_info"),
    ref_info=config.get("ref_info"),
    seq_method=config.get("seq_method"),
    sample_name=config.get("SAMPLE"),
  run:
    sample=list()
    for i,j in zip(range(0,len(SAMPLE),step),GROUP):
      sample+=(f"{j}:",*(f"  - {p}" for p in SAMPLE[i:i+step]))
    sample="\n".join(sample)
    yaml=f"""# common parameter
sample_order:
{group_order}
{sample}
species: {GENOME}
chr_num: 1
fpkm: False
correlation: True 

genome_report_resolution: 1000
bychr_report_resolution: 1000
loop_report_resolution: 1000
TAD_report_method: insulation
TAD_report_resolution: 1000
fire_report_resolution: 1000
compact_report_resolution: 5000
report_name: {params.title}
species_info: {params.species_info}
sample_info: {params.sample_info}
ref_info: {params.ref_info}
project_num: {project}
seq_method: {params.seq_method}
sample_name: {params.sample_name}
compare: {MINUS}
enzyme: MboI
outdir: .
motif: True
TAD_border_annotation: True
target_gene: False
"""
    open(*output,"w").write(yaml)


