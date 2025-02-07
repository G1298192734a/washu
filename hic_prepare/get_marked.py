#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import fire
from bisect import bisect

#pdb文件以N结尾，代表最后一个bin
#开头处有一个N，到N处代表第一个bin
#中间处bin均匀分布

def mark(pdb,fai,bed,res,out_t="bin2pdb.txt",out_r="pdb2bin.txt",out_pdb=False):
    from itertools import chain
    ns,res=list(),int(res) # N位置list 分辨率转数字

    # 分染色体bin总数list,算位置需要
    bsums=[-(-int(i.split()[1])//res) for i in open(fai)]

    # PDB
    # pdb分链list,输出pdb需要
    pdbs=[list() for i in range(len(bsums))]
    for i in open(pdb):
        pdbs[ord(i[-50])-ord("A")].append(i.rstrip())
        # 记录N位置,算位置需要
        if i[-59]=="N":ns.append(int(i.split()[1]))

    bins=[
        [1+int(res/i*p) for p in range(i)]+[1+res+int(res*(k-2)/(j-i-1)*p) for p in range(j-i)]
        for i,j,k in zip(ns[0::2],ns[1::2],bsums)
    ]

    bin_pdbs,pdb_bins=list(),list()
    for i in chain(*pdbs):
        i_l=i.split()
        n,chrom=int(i_l[1]),ord(i_l[4])-ord("A")
        pdb_bins.append([i_l[4],i_l[1],f"chr{chrom+1}",str(bins[chrom][n-1])]+i_l[5:8])
    # bed chrom start end python格式
    for i in open(bed):
        i_l=i[3:].split() #一般chr开头去掉,split会去掉末尾\n
        chrom,start,end=map(int,i_l[0:3])
        chrom-=1 #转坐标
        # 算位置
        s,e=(bisect(bins[chrom],j) for j in (start,end))
        
        # 输出位置映射表(统一用pdb的samtools格式)
        chains=chr(chrom+ord('A'))
        bin_pdbs.append(map(str,[i.rstrip(),chains,s,e]))

        # 如果bed有第四列，把pdb改到对应原子
        if len(i_l)==3:
            continue
        else:
            out_pdb=pdb.replace(".pdb","_marked.pdb")
        atm=i_l[3]
        # 修改pdbs的原子
        for j in range(s-1,e):
            if j<0:continue
            if len(atm)==1: atm=" "+atm
            pdbs[chrom][j]=pdbs[chrom][j].replace(" O",atm)

    # 输出映射表
    with open(out_r,"w") as OUT:
        print(*map("\t".join,pdb_bins),sep="\n",file=OUT)
    with open(out_t,"w") as OUT:
        print(*map("\t".join,bin_pdbs),sep="\n",file=OUT)

    # 输出pdb
    if not out_pdb: return
    with open(out_pdb,"w") as OUT:
        print(*chain(*pdbs),sep="\n",file=OUT)

if __name__=="__main__":
    fire.Fire(mark)


