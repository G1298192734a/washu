#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import fire,re

# 之前我认为从gff获取mapfile.txt会快一些，现在看来是不对劲的
# 现在仅从fa文件生成mapfile.txt
def Map(fa,patt=".",sort=False,outfile="mapfile",ext="mitochondrion"):
  # 匹配pattern
  ids=[i.split()[0][1:] for i in open(fa) if re.match(rf">{patt}",i) and not re.search(rf"{ext}",i)]
  # 排序,必须带有捕获
  if sort:
    ids=sorted(ids,key=lambda x:int(re.search(rf"{sort}",x).group(1) if re.search(rf"{sort}",x) else 100))

  with open(f"{outfile}.txt","w") as OUT:
    print(*(f"{j}\tchr{i+1}" for i,j in enumerate(ids)),sep="\n",file=OUT)
  with open(f"{outfile}_hic.txt","w") as OUT:
    print(*(f"chr{i+1}\tchr{i+1}" for i in range(len(ids))),sep="\n",file=OUT)
  
if __name__=="__main__":
  fire.Fire(Map)


