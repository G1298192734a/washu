import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import fire

def read_csv_chunk(filePath):
  dfs=pd.read_table(filePath,chunksize=1000000,dtype={"0":"str","3":"int","4":"int","5":"str"},usecols=[0,3,4,5],header=None,engine="c",na_filter=False,low_memory=False)
  context_info={i:[] for i in ("CG","CHG","CHH")}
  for p in dfs:
    for i,j in p.groupby(5):
      chr_mc=j.groupby(0)[[3]].sum()
      depth=(j[3]+j[4]).value_counts()
      context_info[i].append((j[3].sum(),depth.sum(),depth[depth.index<=50],chr_mc))

  classifys,chr_clf=[],[]
  for i,j in context_info.items():
    csite=sum((p[1] for p in j))
    classifys.append([i, csite, sum((p[0] for p in j))])
    depth_all=pd.concat([p[2] for p in j],axis=1,join="outer").sum(axis=1).sort_index()
    chr_all=pd.concat([p[3] for p in j],axis=1,join="outer").sum(axis=1).sort_index()
    chr_clf+=[(p,i,int(q)) for p,q in zip(chr_all.index,chr_all)]
    context_info[i]=(csite,depth_all)
  classify_tab=pd.DataFrame(classifys,columns=("ctype","C","mC"))
  classify_tab["mC_rate"]=classify_tab["mC"]/classify_tab["C"]
  classify_tab["ctype_rate"]=classify_tab["mC"]/classify_tab["mC"].sum()
  classify_tab.to_csv("ctype_classify.txt", index=None, sep="\t")

  chr_clf_tab=pd.DataFrame(chr_clf,columns=("chrom","ctype","mC"))
  chr_clf_tab["mC_rate"]=chr_clf_tab["mC"]/pd.concat([chr_clf_tab.groupby("chrom")[["mC"]].sum()]*3)["mC"].tolist()
  chr_clf_tab.to_csv("ctype_classify_chr.txt", index=None, sep="\t")

  return context_info

def plot_coverage(depth):
  plot_data=read_csv_chunk(depth)
  fig, axs = plt.subplots(1, 1)
  for i,j in plot_data.items():
    csite,depth_all=j
    x_smooth=np.linspace(0,50,300)
    def accumu():
      a,b=0,csite
      while True:
        yield b
        a,b=a+1,b-depth_all[a]
    acm=accumu()
    spl=make_interp_spline(depth_all.index,[next(acm)/csite for p in range(51)])
    y_smooth=spl(x_smooth)
    axs.plot(x_smooth,y_smooth,label=i)
  fig.legend(loc=1, bbox_to_anchor=(1,1), bbox_transform=axs.transAxes,fancybox=True, framealpha=0.5)
  plt.savefig("csite_cover_distribution.pdf")

if __name__=="__main__":
  fire.Fire(plot_coverage)


