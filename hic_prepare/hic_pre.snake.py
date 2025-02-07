#/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import yaml
import time

DIR=os.path.abspath(config["DIR"])
#project=f"whfs-xs-{os.path.basename(config['DIR']).split('_',1)[0]}"
project="whfs-xs-232026-01"
#REF=os.path.join(DIR,"00.ref")
#RAW=os.path.join(DIR,"01.rawdata")
#HIC=os.path.join(DIR,"03.HiC")
#DIFF=os.path.join(DIR,"04.HiC_diff")
GENOME=config["GENOME"]
SPECIES=config["SPECIES"]
TITLE=config["TITLE"]
SAMPLES=config["SAMPLES"].split(",")
DATANAMES=config["DATANAMES"].split(",")
MINUS=config["MINUS"]

soft_config = yaml.load(open(config["PRE_SOFT_YAML"], "r"), Loader=yaml.FullLoader)
config = {**config, **soft_config}
python3=config["python3"]
s=" "*4
date=time.strftime("%Y%m%d")

#spls=SAMPLES.split(",")
#dats=DATANAMES.split(",")
#groups=list({j for i in MINUS.split(",") for j in i.split(".minus.")})

rule all:
    input:
#        expand(os.path.join(DIR,project+"_{sample}","00.ref","work.sh"),sample=SAMPLES),
#        expand(os.path.join(DIR,project+"_{sample}","00.ref","{sample}.mrna.fasta"),sample=SAMPLES),
#        expand(os.path.join(DIR,project+"_{sample}","01.hic","single_sample.yaml"),sample=SAMPLES),
#        expand(os.path.join(DIR,project+"_{sample}","01.hic","work.sh"),sample=SAMPLES),
        expand(os.path.join(DIR,project+"_{sample}","01.hic","08.report","report.yaml"),sample=SAMPLES),
#        expand(os.path.join(DIR,project+"_{sample}","01.hic","08.report","work.sh"),sample=SAMPLES),
#        os.path.join(REF,"work.sh"),
#        os.path.join(RAW,"work.sh"),
#        os.path.join(DIR,"02.RNA_result","work.sh"),
#        os.path.join(HIC,"single_sample.yaml"),
#        os.path.join(HIC,"work.sh"),
#        os.path.join(DIFF,"multi_sample.yaml"),
#        os.path.join(DIFF,"work.sh"),
#        os.path.join(REF,GENOME+".mrna.fasta"),
#        os.path.join(RAW,"md5sum.txt"),
#        os.path.join(HIC,"08.report","report.yaml"),
#        os.path.join(HIC,"08.report","work.sh"),
#        os.path.join(DIFF,"05.multi_report","report.yaml"),
#        os.path.join(DIFF,"05.multi_report","work.sh"),

rule ref:
    output:
        os.path.join(DIR,"{project}_{sample}","00.ref","work.sh"),
    params:
        get_map=config["get_map"],
        prepare=config["ref_prepare"],
        sample="{sample}",
        fa="/public/frasergen/backup/lab/2024/whfs-xs-232026/whfs-xs-232026_240301/silkpan/silkworm_genome/{sample}.Final.v1.fa",
        gff="/public/frasergen/backup/lab/2024/whfs-xs-232026/whfs-xs-232026_240301/silkpan/silkworm_genome/{sample}.v1.gff",
        fimo=config.get("fimo"),
    run:
        sh=f"""
#!/usr/bin/bash
# 准备参考和注释
# 生成mapfile
{python3} {params.get_map} gff {params.gff} --pattern Chr
# ref prepare process
module load gffread/0.12.6
export PATH=/public/frasergen/PUB/software/python/Python-3.7.7/bin:$PATH
bash {params.prepare} {params.fa} {params.gff} mapfile.txt {params.fimo} {params.sample}
        """
        with open(*output,"w") as SH:
            print(sh,file=SH)

rule ref_sh:
    output:
        os.path.join(DIR,"{project}_{sample}","00.ref","{sample}.mrna.fasta"),
    params:
        os.path.join(DIR,"{project}_{sample}","00.ref")
    run:
        shell(f"cd {params[0]} && bash work.sh")

#rule rawdata:
#    output:
#        os.path.join(DIR,RAW,"work.sh"),
#    params:
#        rawdata=config["get_rawdata"],
#        sub=config["get_sub"],
#        dt_stg=os.path.join(DIR,config["data_storage"],project,f"{project}_HiC_submit"),
#    run:
#        stg=[os.path.join(DIR,params.dt_stg,i) 
#            for i in ("00.rawdata","01.report","02.result")]
#        snakemake.utils.makedirs(stg)
#        sh=f"""
##!/usr/bin/sh
## 准备数据路径
##{python3} {params.rawdata} 数据源导出.xlsx {stg[0]}
## 拆分脚本
##{python3} {params.sub} --insh cp.sh numsign
### 验证md5
##md5sum -c md5sum.txt > check.txt
### 生成md5
##ln -s {os.path.join(DIR,stg[0],'*gz')} .
##md5sum *gz > md5.txt
##mv md5.txt {stg[0]}
#        """
#        with open(*output,"w") as SH:
#            print(sh,file=SH)
#
#rule cp_sh:
#    output:
#        os.path.join(DIR,RAW,"md5sum.txt"),
##    params:
##        os.path.join(DIR,RAW,"数据源导出.xlsx"),
#    run:
#        shell(f"cd {RAW} && bash work.sh")
#
#rule rna:
#    output:
#        os.path.join(DIR,DIR,"02.RNA_result","work.sh"),
#    run:
#        sh="""
##!/usr/bin/sh
## 准备RNA result <fpkm> <deg>
## fpkm
## DEG
#        """
#        with open(*output,"w") as SH:
#            print(sh,file=SH)

rule single:
    output:
        yaml=os.path.join(DIR,"{project}_{sample}","01.hic","single_sample.yaml"),
        sh=os.path.join(DIR,"{project}_{sample}","01.hic","work.sh"),
    params:
        maps=os.path.join(DIR,"{project}_{sample}","00.ref","mapfile_hic.txt"),
        prefix=os.path.join(DIR,"{project}_{sample}","00.ref","{sample}"),
        hic3="/work/frasergen/3D/work/shaojie/script/HiC/PGT/hic3_workflow.py",
        jaspar=config["fimo"],
        sample="{sample}",
        hic=os.path.join(DIR,"{project}_{sample}","01.hic"),
    run:
        snake=f"snakemake -ps {params.hic3} --configfile single_sample.yaml --cores 10 --cluster 'sbatch -p xhacexclu12,xhacexclu16,xhacexclu03 -n 1 -N 1 --mem 300G --cpus-per-task=1' --jobs 10 --rerun-incomplete"
        sh=f"""
module load bwa/0.7.17 samtools/1.12 python/3.7.10 snakemake perl/5.34.0 circos R/3.6.3 gcc/9.2.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/frasergen/PUB/software/hicexplorer/HiCExplorer-3.7.2/hdf5/lib
# nohup投递
{snake}
        """
        with open(output.sh,"w") as SH:
            print(sh,file=SH)
        yaml=f"""
reference: {params.prefix}.fasta
gene_bed: {params.prefix}_gene.bed
chrom_map: {params.maps}
gene_mrna: {params.prefix}.gene2mark
mrna_fa: {params.prefix}.mrna.fasta
jaspar_db: {params.jaspar}
genome_motif: /work/frasergen/3D/work/shaojie/script/HiC/PGT/whfs-xs-232026-01/whfs-xs-232026-01_{params.sample}/00.ref/fimo.tsv
# species choose from 'animal, plant, fungi, bacteria'
species: a
pattern: Chr
platform: BGI
samplenames: {params.sample}
resolution:
    basic: 5000
    genome: 100000,200000
    compartment: 100000
    tad: 40000
    loop: 10000
    reproduce: 200000
        """
#        rawdatas=[]
#        for i in range(len(spls)):
#            s1,s2=[f"{s}{spls[i]}{j}:" for j in range(1,3)]
#            r1,r2=[os.path.join(DIR,RAW,f"{dats[i]}{j}") for j in range(1,3)]
#            r1=[f"{s*2}- {r1}_R{j}.fq.gz" for j in range(1,3)]
#            r2=[f"{s*2}- {r2}_R{j}.fq.gz" for j in range(1,3)]
#            rawdatas+=[f"{groups[i]}:",s1,*r1,s2,*r2]
        with open(output.yaml,"w") as YAML:
            print(yaml,file=YAML)
        shell(f"cd {params.hic} && {snake} --unlock && nohup bash work.sh &")

#rule diff:
#    output:
#        yaml=os.path.join(DIR,DIFF,"multi_sample.yaml"),
#        sh=os.path.join(DIR,DIFF,"work.sh"),
#    params:
#        hic=os.path.join(DIR,HIC),
#        biorep=config["hic3_biorep"],
#    run:
#        sh=f"""#!/usr/bin/sh
#module load bwa/0.7.17 samtools/1.12 python/3.7.10 snakemake perl/5.34.0 circos R/3.6.3 gcc/9.2.0
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/public/frasergen/PUB/software/hicexplorer/HiCExplorer-3.7.2/hdf5/lib
#export R_LIBS=/public/frasergen/3D/work/wanghao/R_pkgs_add_on/20231006.R-4.2.1_clusterProfiler
#snakemake -ps {params.biorep} --configfile multi_sample.yaml --cores 96 --cluster 'sbatch -p xhacexclu03,xhacexclu12,xhacexclu16 -n 1 -N 1 --cpus-per-task=12' --jobs 8
#        """
#        with open(output.sh,"w") as SH:
#            print(sh,file=SH)
#        yaml=f"""
#single_sample_dir: {params.hic}
## 对照，实验
#sample_order: {','.join(groups)}
## 实验.minus.对照
#sample_minus: {MINUS}
#genome_res: 200000
#compartment_res: 100000
#tad_res: 40000
#loop_res: 10000
#        """
#        with open(output.yaml,"w") as YAML:
#            print(yaml,file=YAML)

rule single_report:
    output:
        yaml=os.path.join(DIR,"{project}_{sample}","01.hic","08.report","report.yaml"),
        sh=os.path.join(DIR,"{project}_{sample}","01.hic","08.report","work.sh"),
    params:
        fai=os.path.join(DIR,"{project}_{sample}","00.ref","{sample}.fasta.fai"),
        sample="{sample}",
        single_pre="/work/frasergen/3D/work/shaojie/script/HiC/PGT/single_report_prepare.py",
        single_report="/work/frasergen/3D/work/shaojie/script/HiC/PGT/single_sample_report",
        hic=os.path.join(DIR,"{project}_{sample}","01.hic"),
        report=os.path.join(DIR,"{project}_{sample}","01.hic","08.report"),
    run:
        g_len=sum((int(i.split()[1]) for i in open(params.fai)))
        yaml=f"""
contract: {project}
title: {TITLE}
species: {params.sample}
reference: {params.sample}
library: Hi-C 2.0建库，内切酶MboI，Illumina PE150测序
chrom_num: 2n=56
genome_size: {g_len:,} bp
compartment: cooltools
tad: cooltools
loop: fithic
basic_res: 5000
genome_res: 100000
reproduce_res: 200000
compartment_res: 100000
tad_res: 40000
loop_res: 10000
report_outdir: prepare_report
samplenames:
    - {params.sample}: {params.sample}
        """
#        sp_names=[f"{s}- {i}: {j}1,{j}2" for i,j in zip(spls,spls)]
        with open(output.yaml,"w") as YAML:
            print(yaml,file=YAML)
        sh=f"""#!/usr/bin/bash
#配置文件必须命名为report.yaml
export PATH=/public/frasergen/PUB/software/python/Python-3.7.7/bin/:$PATH
python3 {params.single_pre} --indir {params.hic} --yaml_file report.yaml
export PATH=/public/frasergen/PUB/software/Anaconda/anaconda3-3d/envs/report/bin/:$PATH
{params.single_report}
zip -r {project}_{params.sample}_single_report_{date}.zip prepare_report single_report.html
zip -r {project}_{params.sample}_single_result_{date}.zip result
        """
        with open(output.sh,"w") as SH:
            print(sh,file=SH)
        shell(f"cd {params.report} && nohup bash work.sh &")

#rule diff_report:
#    output:
#        yaml=os.path.join(DIR,DIFF,"05.multi_report","report.yaml"),
#        sh=os.path.join(DIR,DIFF,"05.multi_report","work.sh"),
#    params:
#        multi_pre=config["multi_pre"],
#        multi_report=config["multi_report"],
#    run:
#        yaml=f"""
#contract: {project}
#title: {TITLE}
#sample_minus: {MINUS}
#genome_res: 200000
#compartment_res: 100000
#tad_res: 40000
#loop_res: 10000
#report_outdir: prepare_report
#multi_sample_dir: {DIFF}
#        """
#        with open(output.yaml,"w") as YAML:
#            print(yaml,file=YAML)
#        sh=f"""#!/usr/bin/bash
##配置文件必须命名为report.yaml
#export PATH=/public/frasergen/PUB/software/python/Python-3.7.7/bin/:$PATH
#python3 {params.multi_pre} --indir {DIFF} --yaml_file report.yaml
#export PATH=/public/frasergen/PUB/software/Anaconda/anaconda3-3d/envs/report/bin/:$PATH
#{params.multi_report}
#zip -r {project}_multi_report_{date}.zip prepare_report multi_sample_report.html
#zip -r {project}_multi_result_{date}.zip result
#        """
#        with open(output.sh,"w") as SH:
#            print(sh,file=SH)


