import yaml
import os
from pathlib import Path
from itertools import chain

DIR=Path(config.get("Dir")).resolve()
dirs=("00.longest_gtf","01.genome","02.single_pre","03.multi_pre")
GTF,GEN,SINGLE,MULTI=map(lambda x:DIR.joinpath(x), dirs)

f=lambda x:yaml.load(open(x,"r"),Loader=yaml.FullLoader).items()
others=dict(chain(*map(f,config.get("YAML"))))
config={**config,**others}

FAI,HiC,Diff=map(lambda x:Path(config.get(x)), ("Fai","HiC","multi_sample_dir"))
prefix=FAI.stem[:-6]
CHR=[i.split()[0] for i in open(FAI)]
washu = config.get("washu_dir")

f=lambda x:config.get(x).split(",")
f1=lambda x:x.keys()
SPL,MINUS=[*chain(*map(f1,config.get("samplenames")))],f("sample_minus")

def get_cool(wildcards):
    return next(HiC.glob(f"06*/01*/{wildcards.spl}*{config.get('tad_res')}*"))
def get_compt(wildcards):
    return next(HiC.glob(f"05*/02*/{wildcards.spl}*{config.get('compartment_res')}/*gc_gene.bed"))
def get_tad(wildcards):
    return next(HiC.glob(f"06*/02*/{wildcards.spl}*{config.get('tad_res')}/TAD_genome.bed"))
def get_loop(wildcards):
    return next(HiC.glob(f"07*/02*/*/{wildcards.spl}.combine.mustache")) if config.get("loop")=="mustache" else                                                     next(HiC.glob(f"07*/03*/{wildcards.spl}*{config.get('loop_res')}/04*/*cis*gz"))
def get_diff(wildcards):
    return next(Diff.glob(f"01*/02*/{wildcards.minus}"))

rule all:
    input:
        GEN.joinpath(washu),
        DIR.joinpath(washu),

rule get_gtf:
    output:
        GTF.joinpath(f"{prefix}_mRNA_genename.gtf"),
    params:
        perl=config.get("perl"),
        longest_gtf=config.get("longest_gtf"),
        gffread=config.get("gffread"),
    run:
        script=rf"""#!/usr/bin/sh
export PERL5LIB=/public/frasergen/3D/work/wanghao/myconda/AGAT/lib/site_perl/5.26.2:$PERL5LIB
{params.perl} {params.longest_gtf} -gff {FAI.parent.joinpath(prefix)}.gff -o {prefix}.gff && \
{params.gffread} {prefix}.gff -T -o {prefix}.gtf && \
awk -F'\t|;| ' -v OFS="\t" 'NF==14{{print $1,$2,$3,$4,$5,$6,$7,$8,$9 " " $13 "; " $12 " " $13 ";";next}}{{print $1,$2,$3,$4,$5,$6,$7,$8,$9 " " $13 "; " $12 " " $13 "; " $15 " " $16 ";"}}' {prefix}.gtf > {output[0]}
        """
        open(GTF.joinpath("work.sh"),"w").write(script)
        shell(f"cd {GTF} && bash work.sh")

rule genome:
    input:
        rules.get_gtf.output,
    output:
        directory(GEN.joinpath(washu)),
    params:
        config.get("genome_prepare"),
    run:
        script=f"""#!/usr/bin/sh
perl {params[0]} -fa {FAI.parent.joinpath(FAI.stem)} -gtf {input[0]} && bash run.sh
/usr/bin/rename genome {{0}} result/* &&  mv result {{0}}
#rsync -P --rsh=ssh -r {{0}} shaojie@114.115.134.103:/local_data1/3D/project/browser/data/config
        """.format(washu)
        open(GEN.joinpath("work.sh"),"w").write(script)
        shell(f"cd {GEN} && bash work.sh")

rule single_pre:
    output:
        directory(SINGLE.joinpath("{spl}")),
        SINGLE.joinpath("{spl}","work.sh"),
    params:
        config.get("bychr"),
        get_cool,
        get_compt,
        get_tad,
        get_loop, 
        "/public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin/pigz",
        "/public/frasergen/PUB/software/bedtools/bedtools-2.30.0/bin/bedtools",
    run:
        x=config.get("loop_res")//2
        loop_sh=f"sed 1d {params[-3]}" if config.get("loop")=="mustache" else                 rf"""{params[-2]} -p 30 -cd {params[-3]}|perl -anle 'print join "\t","chr$F[0]",$F[1]-{x},$F[1]+{x},"chr$F[2]",$F[3]-{x},$F[3]+{x},$F[4]'|{params[-1]} sort -faidx {FAI} -i stdin"""
        script=rf"""#!/usr/bin/sh
mkdir 00.bychr && cd 00.bychr
export PATH=/public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin:$PATH
python {{0}} {{1}} && cd - && mkdir 01.compartment && \
perl -anle 'next if $F[6] ne "A" && $F[6] ne "B";print join "\t",@F[0..3,6]' {{2}} > 01.compartment/sample_100000.compartment_cscore.bedgraph.new && \
mkdir 02.TAD && cp {{3}} 02.TAD && \
mkdir 03.loop && {loop_sh} > 03.loop/cis_interaction.sort
        """.format(*params[:-3])
        open(output[1],"w").write(script)
        shell(f"cd {output[0]} && bash work.sh")

rule single:
    input:
        rules.single_pre.output,
    output:
        directory(DIR.joinpath("{spl}")),
        DIR.joinpath("{spl}","work.sh"),
    wildcard_constraints:
         spl="|".join(SPL),
    params:
        config.get("single_prepare"),
    run:
        bychr,compt,tad,loop=map(lambda x:Path(input[0]).joinpath(*x), 
        (("00.bychr",),("01.compartment","sample_100000.compartment_cscore.bedgraph.new"),
              ("02.TAD","TAD_genome.bed"),("03.loop","cis_interaction.sort"))
        )
        script=f"""#!/usr/bin/sh
perl {params[0]} -fai {FAI} -hm_path {bychr} -hm_res {config.get("tad_res")} -compt_path {compt} -compt_res {config.get("compartment_res")} -tad_path {tad} -tad_res {config.get("tad_res")} -loop_path {loop} -s_id {wildcards.spl} && bash run.sh
        """
        open(output[1],"w").write(script)
        shell(f"cd {output[0]} && bash work.sh")

rule multi_pre:
    output:
        directory(MULTI.joinpath("{minus}")),
        MULTI.joinpath("{minus}","work.sh"),
    params:
        get_diff,
    run:
        f=lambda x: rf"""
zcat {params[0]}/{wildcards.minus}.{x}.matrix.gz | awk -F'\t' 'NR>1{{for(i=2;i<NF;i++){{printf("%.8f ",$i==""?0:$i)}};printf("%.8f\n",$NF==""?0:$NF)}}' > {x}.minus_heatmap.matrix
        """
        with open(output[1],"w") as OUT:
            print(*map(f,CHR),sep="",file=OUT)
        shell(f"cd {output[0]} && bash work.sh")
#        shell(f"cd {output[0]} && /public/frasergen/PUB/software/python/Python-3.7.10/bin/python3 /work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/get_subscript.py --insh work.sh nrows 1")

rule multi:
    input:
        rules.multi_pre.output,
    output:
        directory(DIR.joinpath("{minus}")),
        DIR.joinpath("{minus}","work.sh"),
    wildcard_constraints:
         minus="|".join(MINUS),
    params:
        config["multi_prepare"],
    run:
        bychr=MULTI.joinpath(wildcards.minus)
        script=f"""#!/usr/bin/sh
perl {params[0]} -fai {FAI} -minus_path {bychr} -minus_res {config.get("genome_res")} -s_id {wildcards.minus} && bash run.sh
        """
        open(output[1],"w").write(script)
        shell(f"cd {output[0]} && bash work.sh")

rule transform:
    input:
        expand(DIR.joinpath("{spl}"),spl=SPL),
        expand(DIR.joinpath("{minus}"),minus=MINUS),
    output:
        directory(DIR.joinpath(washu)),
    run:
       dirs=SPL+MINUS
       script=f"""#!/usr/bin/sh
mkdir {{0}} && cp -r {" ".join(dirs)} {{0}}
rm -r $(ls {{0}}/*/* -d|perl -ne 'print if !/[.*gz|.*tbi]$/')
#rsync -P --rsh=ssh -r {{0}} shaojie@114.115.134.103:/local_data1/3D/project/browser/data/customerdata
       """.format(washu)
       open(DIR.joinpath("work.sh"),"w").write(script)
       shell(f"cd {DIR} && bash work.sh")


