import os,yaml,shutil
#import ffmpeg
from pathlib import Path
from itertools import chain

Sample,Res=map(config.get, ("Heatmap","Res"))
Fai=Path(config["Fai"])
Dir=Path(config.get("Dir")).resolve()
sub_dir=("01.matrix_fill","02.3d_model","03.gif")
Matrix,Model,Gif=map(Dir.joinpath,sub_dir)

soft_config = yaml.load(open(config["Software_yaml"], "r"), Loader=yaml.FullLoader)
config = {**config, **soft_config}

chains=list(map(chr, (*range(65,91),*range(97,123))))
chroms=[i.split()[0] for i in open(Fai)]
rgbs=[i.split()[-1] for i in open(config.get("pml")) if i.startswith("set_color")]

def get_hm(wildcards):
    return Sample[wildcards.sample][0]

rule all:
    input:
        expand(Gif.joinpath("{sample}","chrs","{chrom}","{chrom}_2.png"),sample=Sample.keys(),chrom=chroms),
        expand(Gif.joinpath("{sample}","{sample}.gif"),sample=Sample.keys()),
        Dir.joinpath(re.sub("^analysis",f"whfs-xs-{Fai.parents[2].name[:6]}",config["Dir"]))

rule fill_pre:
    output:
        Matrix.joinpath("bins.bed"),
    run:
        script=f"""#!/usr/bin/sh
bedtools makewindows -g {Fai} -w {Res} > {output[0]}
        """
        open(Matrix.joinpath("work.sh"),"w").write(script)
        shell(f"bash {Matrix.joinpath('work.sh')}")

rule matrix:
    output:
        Matrix.joinpath("{sample}","{sample}_filled_matrix.txt"),
    params:
        get_hm,
        config.get("hm2matrix"),
        config.get("cool2matrix"),
        config.get("matrix_fill"),
    run:
        subdir=Matrix.joinpath(wildcards.sample)
        script=f"""#!/usr/bin/sh
cp {params[0]} . && python3 {params[2] if params[0].endswith("cool") else params[1]} {os.path.basename(params[0])} && cat *txt|awk -v res={Res} -f {params[3]} {Fai} - > {output[0]}
        """
        open(subdir.joinpath("work.sh"),"w").write(script)
        shell(f"cd {subdir} && bash work.sh")

rule pastis:
    input:
        rules.matrix.output,
        rules.fill_pre.output,
    output:
        Model.joinpath("{sample}","{sample}.pse"),
        Model.joinpath("{sample}","{sample}.pml"),
    params:
        config.get("python_pastis").get("path"),
        config.get("hm2heat"),
        config.get("legend"),
        "/work/frasergen/3D/work/shaojie/script/HiC/3d_model/get_marked.py",
    run:
        subdir=Model.joinpath(wildcards.sample)
        script=f"""#!/usr/bin/sh
{params[0]} {params[1]} -i {Fai} -m {input[0]} -r {Res} -s {Res//1000} && cp *{Res}.pdb {wildcards.sample}.pdb

cat <<END > {wildcards.sample}.pml
load {wildcards.sample}.pdb
bg white
show_as sphere
orient
set ray_shadows, 0
set sphere_scale, 0.6
{"{0}".join((f"select {i}, chain {j}" for i,j in zip(chroms,chains)))}
{"{1}".join((f"set_color c{i}, {j}" for i,j in zip(chroms,rgbs)))}
{"{2}".join((f"color c{i}, {i}" for i in chroms))}
select ter, name N
set sphere_scale, 1.0, ter
orient ter
zoom all,2,complete=1
set_color cter, [0,0,255]
color cter, ter
save {output[0]}
END

python3 {params[2]} {wildcards.sample}.pml
python3 {params[3]} {wildcards.sample}.pdb {Fai} {input[1]} {Res}
        """.format("\n","\n","\n")
        open(subdir.joinpath("work.sh"),"w").write(script)
        shell(f"cd {subdir} && bash work.sh")

rule chrom:
    input:
        Model.joinpath("{sample}","{sample}.pse"),
        Model.joinpath("{sample}","{sample}.pml"),
    output:
        Gif.joinpath("{sample}","chrs","{chrom}","{chrom}_2.png"),
    params:
        config.get("pymol").get("path"),
    run:
        subdir=Gif.joinpath(wildcards.sample,"chrs",wildcards.chrom)
        script=f"""#!/usr/bin/sh
cat <<END > {wildcards.chrom}.pml
load {input[0]}
select ter, {wildcards.chrom} & name N
hide
show sphere, {wildcards.chrom}
orient ter
rotate z,90,{wildcards.chrom}
zoom {wildcards.chrom},5,complete=1
png {wildcards.chrom}_1,1000,1000,300,1,1
rotate x,180,{wildcards.chrom}
png {wildcards.chrom}_2,1000,1000,300,1,1
END

{params[0]} -c {wildcards.chrom}.pml
        """
        open(subdir.joinpath("work.sh"),"w").write(script)
        shell(f"cd {subdir} && bash work.sh")

rule png:
    input:
        Model.joinpath("{sample}","{sample}.pse"),
        Model.joinpath("{sample}","{sample}.pml"),
    output:
        Gif.joinpath("{sample}","pngs","fig_{batch}.png"),
    params:
        config.get("pymol").get("path"),
    run:
        subdir=Gif.joinpath(wildcards.sample,"pngs")
        batch=int(wildcards.batch)
        dim,angle="x" if batch<=120 else "y",(batch-10)%120*3
        pngs=(f"png fig_{i:03},1000,1000,300,1,1" for i in range(batch-9,batch+1))
        batch//=10
        script=f"""#!/usr/bin/sh
cat <<END > gif_{batch:03}.pml
load {input[0]}
rotate {dim}, {angle}
{"{0}".join(chain(*([i,j] for i,j in zip([f"rotate {dim}, 3"]*10,pngs))))}
END

{params[0]} -c gif_{batch:03}.pml
        """.format("\n")
        open(subdir.joinpath(f"work_{batch:03}.sh"),"w").write(script)
        shell(f"cd {subdir} && bash work_{batch:03}.sh")

Batch=[f"{i*10:03}" for i in range(1,25)]
rule gif:
    input:
        expand(Gif.joinpath("{{sample}}","pngs","fig_{batch}.png"),batch=Batch),
    output:
        Gif.joinpath("{sample}","{sample}.gif"),
        Gif.joinpath("{sample}","{sample}.mp4"),
    params:
        config.get("gif"),
        config.get("ffmpeg"),
    run:
        subdir=Gif.joinpath(wildcards.sample)
        with open(subdir.joinpath("pngs","png.txt"),"w") as PNG:
            print(*(f"pngs/fig_{i:03}.png" for i in range(1,241)),sep="\n",file=PNG)
#        ffmpeg.input(subdir.joinpath("pngs","*.png"),pattern_type='glob',framerate=25).output(output[1],c="libx264",pix_fmt="yuv420p",r=29.97,crf=18).run()
        script=f"""#!/usr/bin/sh
python3 {params[0]} {subdir.joinpath("pngs","png.txt")} {output[0]}
{params[1]} -pattern_type glob -i "{subdir.joinpath("pngs","*.png")}" -c:v libx264 -pix_fmt yuv420p -crf 18 -r 29.97 {output[1]}
        """
        open(subdir.joinpath("work.sh"),"w").write(script)
        shell(f"cd {subdir} && bash work.sh")

rule result:
    input:
        expand(rules.gif.output,sample=Sample.keys())
    output:
        directory(Dir.joinpath(re.sub("^analysis",f"whfs-xs-{Fai.parents[2].name[:6]}",config["Dir"])))
    run:
      shutil.copytree(Gif,output[0])
      for i in [*Model.glob("*/legend*")][:2]: shutil.copy(i,output[0])
      for i in Path(output[0]).rglob("*"):
          if not re.search("pml|sh$|txt$",i.name): continue
          if i.is_dir(): i.rmdir()
          else: i.unlink()
      shell("find {output} -name .snakemake_timestamp|xargs -I [] rm []")



