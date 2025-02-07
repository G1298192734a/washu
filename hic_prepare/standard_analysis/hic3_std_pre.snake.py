#/usr/bin/env python3
# -*- coding: utf-8 -*-
import yaml,time,re
from pathlib import Path
from itertools import chain

DIR=Path(config["DIR"]).resolve()
#project="whfs-xs-{0}".format(re.match(r"\d+",DIR.name)[0])
project="whfs-xs-{0}".format(re.search(r"(?<=/)\d{6}(?=_)",str(DIR))[0])
dirnames=("00.ref","01.rawdata","02.personal","03.HiC","04.HiC_diff",)
REF,RAW,RNA,HIC,DIFF=map(DIR.joinpath,dirnames)
SAMPLE,DATANAME=map(lambda x:config.get(x).split(","),("SAMPLE","DATANAME"))
GENOME,kingdom,platform=map(config.get,("GENOME","kingdom","library"))
MINUS=config.get("MINUS").split(",")
group=[j for i in MINUS for j in i.split(".minus.")]
GROUP=sorted(set(group),key=group.index) if group else SAMPLE
step=len(SAMPLE)//len(GROUP)

dirs=("01.genome_gene_express","02.diff_gene_express","03.BAM_2_coverage")
fpkm,deg,bam=map(lambda x:RNA.joinpath(x),dirs)

jaspar_dict={
  "a":"JASPAR2022_CORE_vertebrates_non-redundant_v2.meme",
  "p":"JASPAR2022_CORE_plants_non-redundant_v2.meme",
  "f":"JASPAR2022_CORE_fungi_non-redundant_v2.meme",
}

soft_config = yaml.load(open(config["PRE_SOFT_YAML"], "r"), Loader=yaml.FullLoader)
config = {**config, **soft_config}
date=time.strftime("%Y%m%d")

rule all:
  input:
    *(i.joinpath("work.sh") for i in (REF,RAW,RNA,HIC,DIFF)),
    HIC.joinpath("08.report","work.sh"),
    DIFF.joinpath("05.multi_report","work.sh"),

rule ref:
  output:
    REF.joinpath("work.sh"),
  params:
    "/public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/parallel_fimo.py",
    get_map=config["get_map"],
    prepare=config["ref_prepare"],
    jaspar=Path(config["jaspar_db"]).joinpath(jaspar_dict[kingdom]),
  run:
    script=f"""#!/usr/bin/bash

# 准备参考和注释
 fa= #raw fasta
gff= #raw gff
  p= #chrom pattern
  s= #sort pattern
  j={params.jaspar} # fimo

# 生成mapfile
export PATH=/public/frasergen/PUB/software/python/Python-3.7.7/bin:$PATH
python3 {params.get_map} $fa "$p" "$s"

# ref prepare process
module load gffread/0.12.6
bash {params.prepare} $fa $gff mapfile.txt {GENOME}
#python3 {params[0]} --infasta {GENOME}.fasta --meme_db $j --outdir ./
    """
    open(*output,"w").write(script)

rule rawdata:
  output:
    RAW.joinpath("work.sh")
  params:
    get_rawdata=config["get_rawdata"],
    get_sub=config["get_sub"],
    data_storage=config.get("data_storage"),
  run:
    main=f"{project}_HiC_submit_{date}"
    sub="{00.rawdata,01.report,02.result}"
    f=lambda x:Path(params.data_storage).joinpath(project,main,*x)
    data_dir,submit_dir,data=map(f,(
      ["00.rawdata"], [sub], ["00*","*gz"],
    ))
    script=f"""#!/usr/bin/sh

# 准备数据路径
x=数据源导出.xlsx #默认从系统中下载数据文档
mkdir -p {submit_dir}
export PATH=/public/frasergen/PUB/software/python/Python-3.7.10/bin:$PATH
python3 {params.get_rawdata} $x

# 拆分脚本
#python3 {params.get_sub} --insh cp.sh numsign

## 生成备份数据 md5
#cat *md5.txt > md5.txt && rm *_md5.txt
#mv md5.txt *gz {data_dir} && ln -s {data} .

## 验证原始数据 md5,需要在合并完成后再次投递
#md5sum -c md5sum.txt > check.txt
    """
    open(*output,"w").write(script)

rule rna:
  log:
    RNA.joinpath("work.sh")
  params:
    pgt_snake="/work/frasergen/3D/work/shaojie/script/HiC/integrate/pgt_snake.py",
    washu_soft="/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/upload_software.yaml",
    washu_snake="/work/frasergen/3D/work/shaojie/script/HiC/washU/upload/upload_hic3.0_snake.py",
    trid_soft="/work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_soft.yaml",
    trid_snake="/work/frasergen/3D/work/shaojie/script/HiC/3d_model/trid_model_snake.py",
    fpkm_pre=config["fpkm_pre"],
    fpkm=config["fpkm"],
    logfc=config["logfc"],
    integrate="/work/frasergen/3D/work/shaojie/script/HiC/integrate/integrate_snake.py",
  run:
    yamls="\n".join((f"    - {i}" for i in (params.washu_soft,HIC.joinpath("08.report","report.yaml"),DIFF.joinpath("05.multi_report","report.yaml"))))
    cools="\n".join(chain(*((f"    {i}:","    - "+str(HIC.joinpath("05.compartment","01.prepare",f"{i}_100000.cool"))) for i in GROUP)))
    script=f"""#!/usr/bin/sh
#module load python/3.7.10 snakemake
<< 多组学
cat <<END >integrate.yaml
Dir: analysis_HiC_ATAC_RNA
Gene_bed: {HIC.joinpath("01.ref","genome_gene.bed")}
Fpkm: whfs*/02*/gene.fpkm_matrix*
DEG: whfs*/06*/genes/*known.DEG*
Filter_before: whfs*/06*/genes/*before-filter.xls
Spl_vs: rna_cont-vs-treat 
Spls:
{{0}}
ATAC:
  HiC_ATAC:
  - /work/frasergen/backup/3d/project/ATAC/2024
  - atac_cont-vs-treat
YAML: {DIFF.joinpath("05.multi_report","report.yaml")}
END
snakemake -ps {params.integrate} --configfile integrate.yaml -j 10 --cluster 'sbatch -p xhacexclu03,xhacexclu12,xhacexclu16 -n 1 -c 1' --jobs 10
多组学

<< 三维模型
cat <<END >trid.yaml
Fai: {HIC.joinpath("01.ref","genome.fa.fai")}
Dir: analysis_3d_model
Res: 100000
Heatmap:
{cools}
Software_yaml: {params.trid_soft}
END
snakemake -ps {params.trid_snake} --configfile trid.yaml --cores 24 --cluster "sbatch -p xhacexclu16,xhacexclu03,xhacexclu12 -n 1 -c 1" --jobs 24 --rerun-incomplete
三维模型

<< washU
cat <<END >washu.yaml
# 同目录下有gff
Fai: {REF.joinpath(GENOME+".fasta.fai")}
Dir: analysis_washu
washu_dir: {re.match(r".+(?=_)",DIR.name)[0]}
HiC: {HIC}
YAML:
{yamls}
END
snakemake -ps {params.washu_snake} --configfile washu.yaml --cores 24 --cluster "sbatch -p xhacexclu12,xhacexclu16,xhacexclu03 -n 1 -c 1' --jobs 24 --rerun-incomplete
washU

<< PGT
cat <<END >pgt.yaml
Dir: analysis_PGT
HiC: {HIC}
# 同目录下有gff
Fai: {REF.joinpath(GENOME+".fasta.fai")}
YAML: {DIFF.joinpath("05.multi_report","report.yaml")}
BW:
END
snakemake -ps {params.pgt_snake} --configfile pgt.yaml --cores 24 --cluster "sbatch -p xhacexclu12,xhacexclu16,xhacexclu03 -n 1 -c 1" --jobs 24 --rerun-incomplete
PGT

<< hic文件
###python3 /work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/get_subscript.py --insh hic.sh nrows 1 -p xhacexclu12
###/public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin/python /public/frasergen/3D/work/wanghao/my_scripts/python/hiclib_frag_2_pairs/hiclib_frag_2_bgzip_pairs_v0.1.py genome-all-MboI_refined.frag names.pairs genomePath
mkdir java_tmp analysis_hic_file
for names in {" ".join(GROUP)};do 
cat <<END >>hic.sh
/usr/bin/java -Djava.awt.headless=true -Djava.io.tmpdir=java_tmp -Xmx360g -jar /public/frasergen/3D/work/wanghao/soft-sources/github_sources/Juicer_release_1.6/juicer-1.6/scripts/common/juicer_tools.jar pre -j 4 {HIC.joinpath("03.align",'$names".merged.pair.gz"')} analysis_hic_file/$names".hic" {HIC.joinpath("01.ref","genome.chrsize")}
END
done
hic文件
    """.format("\n".join((f"  {i}: {i}" for i in GROUP)))
    open(log[0],"w").write(script)

rule single:
  output:
    yaml=HIC.joinpath("single_sample.yaml"),
    sh=HIC.joinpath("work.sh"),
  params:
    maps=REF.joinpath("mapfile_hic.txt"),
    prefix=REF.joinpath(GENOME),
    hic3=config["hic3"],
    jaspar=Path(config["jaspar_db"]).joinpath(jaspar_dict[kingdom]),
    fimo=REF.joinpath("fimo.tsv"),
  run:
    script=f"""#!/usr/bin/sh
module load bwa/0.7.17 samtools/1.12 python/3.7.10 snakemake perl/5.34.0 circos R/3.6.3 gcc/9.2.0 bedtools
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/frasergen/PUB/software/hicexplorer/HiCExplorer-3.7.2/hdf5/lib
#python3 /work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/get_kegg.py {GENOME} --species {kingdom}
# nohup投递
snakemake -ps {params.hic3} --configfile single_sample.yaml --cores 32 --cluster 'sbatch -p xhacexclu12,xhacexclu16,xhacexclu03 -n 1 -c 8' --jobs 4 --until aln_pair
    """
    open(output.sh,"w").write(script)
    yaml=f"""reference: {params.prefix}.fasta
gene_bed: {params.prefix}_gene.bed
chrom_map: {params.maps}
gene_mrna: {params.prefix}.gene2mark
mrna_fa: {params.prefix}.mrna.fasta
jaspar_db: {params.jaspar}
genome_motif: {params.fimo}
# species choose from 'animal, plant, fungi, bacteria'
species: {kingdom}
platform: {platform}
samplenames: {','.join(GROUP)}
resolution:
  basic: 5000
  genome: 100000,200000
  compartment: 100000
  tad: 40000
  loop: 5000,10000,20000
  reproduce: 200000
    """
    rawdatas=list()
    f=lambda x:(f"- {RAW.joinpath(x)}_R{i}.fq.gz" for i in (1,2))
    for i,j in zip(range(0,len(SAMPLE),step),GROUP):
      a=((f"{p}:",*f(q)) for p,q in zip(SAMPLE[i:i+step],DATANAME[i:i+step]))
      rawdatas+=(f"{j}:",*(f"  {p}" for p in chain(*a)))

    with open(output.yaml,"w") as YAML:
      print(yaml,*rawdatas,sep="\n",file=YAML)

rule diff:
  output:
    yaml=DIFF.joinpath("multi_sample.yaml"),
    sh=DIFF.joinpath("work.sh"),
  params:
    biorep=config["hic3_biorep"],
    biorep1=config["hic3_biorep1"],
  run:
    script=f"""#!/usr/bin/sh
module load bwa/0.7.17 samtools/1.12 python/3.7.10 snakemake perl/5.34.0 circos R/3.6.3 gcc/9.2.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/frasergen/PUB/software/hicexplorer/HiCExplorer-3.7.2/hdf5/lib
export R_LIBS=/public/frasergen/3D/work/wanghao/R_pkgs_add_on/20231006.R-4.2.1_clusterProfiler
snakemake -ps {params.biorep1} --configfile multi_sample.yaml --cores 8 --cluster 'sbatch -p xhacexclu03,xhacexclu16 -n 1 -c 8' --jobs 1 --until diff_matrix
#snakemake -ps {params.biorep1.replace("biorep","biorep1")} --configfile multi_sample.yaml --cores 12 --cluster 'sbatch -p xhacexclu03,xhacexclu12,xhacexclu16 -n 1 -N 1 --cpus-per-task=1' --jobs 12
#snakemake -ps {params.biorep} --configfile multi_sample.yaml --cores 12 --cluster 'sbatch -p xhacexclu03,xhacexclu12,xhacexclu16 -n 1 -N 1 --cpus-per-task=1' --jobs 12
    """
    open(output.sh,"w").write(script)
    yaml=f"""single_sample_dir: {HIC}
# 对照，实验
sample_order: {','.join(GROUP)}
# 实验.minus.对照
sample_minus: {','.join(MINUS)}
genome_res: 200000
compartment_res: 100000
tad_res: 40000
loop_res: 5000
    """
    open(output.yaml,"w").write(yaml)

rule single_report:
  output:
    yaml=HIC.joinpath("08.report","report.yaml"),
    sh=HIC.joinpath("08.report","work.sh"),
  params:
    single_report=config["single_report"],
    single_pre="/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/single_report_prepare.py",
    title=config["TITLE"],
    species=config["SPECIES"],
    loop=config["LOOP"],
  run:
    yaml=f"""contract: {project}
title: {params.title}
species: {params.species}
reference: {GENOME}
library: Hi-C 2.0建库，内切酶MboI，{platform} PE150测序
chrom_num: 2n=
genome_size:  bp
compartment: cooltools
tad: cooltools
loop: {params.loop}
basic_res: 5000
genome_res: 100000
reproduce_res: 200000
compartment_res: 100000
tad_res: 40000
loop_res: 5000
report_outdir: prepare_report
    """
    samples=list()
    for i,j in zip(range(0,len(SAMPLE),step),GROUP):
      samples.append(f"  - {j}: "+",".join(SAMPLE[i:i+step]))
    with open(output.yaml,"w") as YAML:
      print(yaml,"samplenames:",*samples,sep="\n",file=YAML)
    script=f"""#!/usr/bin/bash

#配置文件必须命名为report.yaml
export PATH=/public/frasergen/PUB/software/python/Python-3.7.7/bin/:$PATH
python3 {params.single_pre} --indir {HIC} --yaml_file report.yaml
export PATH=/public/frasergen/PUB/software/Anaconda/anaconda3-3d/envs/report/bin/:$PATH
{params.single_report}
zip -r {project}_single_report_{date}.zip prepare_report single_sample_report.html
zip -r {project}_single_result_{date}.zip result
    """
    open(output.sh,"w").write(script)

rule diff_report:
  output:
    yaml=DIFF.joinpath("05.multi_report","report.yaml"),
    sh=DIFF.joinpath("05.multi_report","work.sh"),
  params:
    multi_pre="/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/multi_report_prepare.py",
    multi_report=config["multi_report"],
    title=config["TITLE"],
    loop=config["LOOP"],
  run:
    yaml=f"""contract: {project}
title: {params.title}
sample_minus: {",".join(MINUS)}
genome_res: 200000
compartment_res: 100000
tad_res: 40000
loop_res: 5000
loop: {params.loop}
report_outdir: prepare_report
multi_sample_dir: {DIFF}
    """
    open(output.yaml,"w").write(yaml)
    script=f"""#!/usr/bin/bash

#配置文件必须命名为report.yaml
export PATH=/public/frasergen/PUB/software/python/Python-3.7.7/bin/:$PATH
python3 {params.multi_pre} --indir {DIFF} --yaml_file report.yaml
export PATH=/public/frasergen/PUB/software/Anaconda/anaconda3-3d/envs/report/bin/:$PATH
{params.multi_report}
zip -r {project}_multi_report_{date}.zip prepare_report multi_sample_report.html
zip -r {project}_multi_result_{date}.zip result
    """
    open(output.sh,"w").write(script)


