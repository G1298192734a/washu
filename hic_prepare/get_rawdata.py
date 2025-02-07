#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import fire,re
from pandas import read_excel
from pathlib import Path
from itertools import chain

def get_raw(spl_x,spl_t="sample",md5_t="md5sum",sh="cp",outdir="./"):
  df=read_excel(spl_x,usecols=["样本名","数据路径"])
  md5s,spl2fq = list(),dict()
  f=lambda x:re.sub(f"[-_]","",x)
  for i in df.itertuples():
    spl,path=i[1:]
    files = [*map(str,Path(path).rglob("*"))]
    spl2fq[spl]=spl2fq.get(spl,[])+[i for i in files if i.endswith("gz") and f(spl) in f(i)]
    md5=chain(*(open(i).read().splitlines() for i in files if "md5" in i.lower() and re.search(r"txt$",i)))
    md5s+=[re.sub(r"\s.*[\s/]",f" {spl}/",i) for i in md5 if f(spl) in f(i)]
  
  with open(f"{md5_t}.txt","w") as MD5:
    print(*sorted(set(md5s),key=lambda x:x.split()[1]),sep="\n",file=MD5)

  with open(f"{sh}.sh","w") as SH,open(f"{spl_t}.txt","w") as SPL:
    for i,j in sorted(spl2fq.items(), key=lambda x:x[0]):
      path=sorted(set(j))
      print(i,*path,sep="\n",file=SPL)
      prefix=Path(outdir).joinpath(i)
      path1,path2=path[0::2],path[1::2]
      print(f"# {i}",f"mkdir {i}",sep="\n",file=SH)
      print("ln","-s",*path1,*path2,i,file=SH)
      print("cat",*path1,">",f"{prefix}_R1.fq.gz",file=SH)
      print("cat",*path2,">",f"{prefix}_R2.fq.gz",file=SH)
      print("md5sum",f"{prefix}*gz",">",f"{prefix}_md5.txt",file=SH)

if __name__=="__main__":
  fire.Fire(get_raw)


