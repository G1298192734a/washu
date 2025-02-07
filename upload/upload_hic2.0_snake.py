import yaml
import os
from pathlib import Path
from itertools import chain

DIR=Path(config.get("DIR")).resolve()
dirs=("00.longest_gtf","01.genome","02.single_pre","03.multi_pre")
GTF,GEN,SINGLE,MULTI=map(lambda x:DIR.joinpath(x), dirs)

FAI,HiC=Path(config.get("FAI")),Path(config.get("HiC"))
prefix=FAI.stem[:-3]
washu=config.get("washu_dir")

f=lambda x:yaml.load(open(x,"r"),Loader=yaml.FullLoader).items()
others=dict(chain(*map(f,config.get("YAML"))))
config={**config,**others}

f=lambda x:config.get(x).split(",")
f1=lambda x:config.get(x).keys()
SPL,MINUS=([*chain(*map(f1,config.get("sample_order")))],f("samplevs"))

def get_heat(wildcards):
    return next(HiC.glob(f"{wildcards.spl}/04*/*/by*/*{config.get('bychr_resolution')//1000}*/heat*"))
def get_compt(wildcards):
    return next(HiC.glob(f"{wildcards.spl}/06*/02*/*{config.get('compartment_resolution')//1000}*/*new"))
def get_tad(wildcards):
    return next(HiC.glob(f"{wildcards.spl}/07*/02*/*{config.get('TAD_resolution')//1000}*/*genome*"))
def get_loop(wildcards):
    return next(HiC.glob(f"{wildcards.spl}/05*/*{config.get('loop_resolution')//1000}*/02*/c*sort"))
def get_diff(wildcards):
    return next(HiC.glob(f"diff/out*/{wildcards.minus}/01*/by*/{wildcards.minus}"))

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
{params.perl} {params.longest_gtf} -gff {FAI.parent.joinpath(FAI.stem[:-3])}.gff -o {prefix}.gff && \
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
/usr/bin/rename {prefix} {{0}} result/* &&  mv result {{0}}
#scp -r {{0}} shaojie@119.3.200.180:/local_data1/3D/project/browser/data/config
        """.format(washu)
        open(GEN.joinpath("work.sh"),"w").write(script)
        shell(f"cd {GEN} && bash work.sh")

rule single:
    output:
        directory(DIR.joinpath("{spl}")),
        DIR.joinpath("{spl}","work.sh"),
    wildcard_constraints:
         spl="|".join(SPL),
    params:
        config.get("genome_prepare"),
        get_heat,
        config.get("bychr_resolution"),
        get_compt,
        config.get("compartment_resolution"),
        get_tad,
        config.get("TAD_resolution"),
        get_loop,
        config.get("loop_resolution"),
    run:
        script=f"""#!/usr/bin/sh
perl {{0}} -fai {FAI} -hm_path {{1}} -hm_res {{2}} -compt_path {{3}} -compt_res {{4}} -tad_path {{5}} -tad_res {{6}} -loop_path {{7}} -loop_res {{8}} -s_id {wildcards.spl} && bash run.sh
        """.format(*params)
        open(output[1],"w").write(script)
        shell(f"cd {output[0]} && bash work.sh")

rule multi:
    output:
        directory(DIR.joinpath("{minus}")),
        DIR.joinpath("{minus}","work.sh"),
    wildcard_constraints:
        minus="|".join(MINUS),
    params:
        config.get("genome_prepare"),
        get_diff,
        config.get("bychr_resolution"),
    run:
        script=f"""#!/usr/bin/sh
perl {params[0]} -fai {FAI} -minus_path {params[1]} -minus_res {params[2]} -s_id {wildcards.minus} && bash run.sh
        """
        open(output[1],"w").write(script)
        shell(f"cd {output[0]} && bash work.sh")

rule transform:
    input:
        (expand(DIR.joinpath("{spl}"),spl=SPL),expand(DIR.joinpath("{minus}"),minus=MINUS),)
    output:
        directory(DIR.joinpath(washu)),
    run:
       dirs=SPL+MINUS
       script=f"""#!/usr/bin/sh
mkdir {{0}} && cp -r {" ".join(dirs)} {{0}}
rm -r $(ls {{0}}/*/* -d|perl -ne 'print if !/[.*gz|.*tbi]$/')
#scp -r {{0}} shaojie@119.3.200.180:/local_data1/3D/project/browser/data/customerdata
       """.format(washu)
       open(DIR.joinpath("work.sh"),"w").write(script)
       shell(f"cd {DIR} && bash work.sh")


