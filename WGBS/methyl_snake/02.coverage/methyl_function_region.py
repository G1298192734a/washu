from pandas import read_table
from numpy import linspace
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import make_interp_spline
import fire,re

def plot_methyl(matrix):
  matrix=read_table(matrix,comment="@",dtype={"3":"str",**{i:"float" for i in range(6,66)}},usecols=(3,*range(6,66)),header=None,engine="c",na_filter=False,low_memory=False)
  matrix[3]=matrix[3].apply(lambda x:re.sub("\d+$","",x))
  plot_data={i:[] for i in ("CG","CHG","CHH")}
  labels=[]
  for i,j in matrix.groupby(3,sort=False):
    labels.append(i)
    mlevel=j.loc[:,6:65].mean()
    plot_data["CG"]+=mlevel.loc[6:25].tolist()
    plot_data["CHG"]+=mlevel.loc[26:45].tolist()
    plot_data["CHH"]+=mlevel.loc[46:65].tolist()

  fig, axs = plt.subplots(3, 1)
  n=0
  for i,j in plot_data.items():
    x_smooth=linspace(0,20*len(labels),120*len(labels))
    spl=make_interp_spline(range(0,20*len(labels)),j)
    y_smooth=spl(x_smooth)
    axs[n].plot(x_smooth,y_smooth)
    axs[n].set_xticks([])
    axs[n].set_title(i,fontdict={"fontsize":10})
    axs[n].set_ylim(axs[n].get_ylim()[0],axs[n].get_ylim()[1])
    axs[n].set_xlim(0,20*len(labels))
    yformatter=ScalarFormatter(useMathText=True)
    yformatter.set_powerlimits((-2,2))
    axs[n].yaxis.set_major_formatter(yformatter)
    for p,q in enumerate(labels):
      if p>0:
        axs[n].vlines(p*20,axs[n].get_ylim()[0],axs[n].get_ylim()[1],color="r" if p==1 else "b" if p==5 else "grey",linestyles="dashdot")
      if n==2:
        axs[n].vlines(10+p*20,0,-0.04,color="black",transform=axs[n].get_xaxis_transform(),clip_on=False)
        axs[n].text(10+p*20,axs[n].get_ylim()[0], q, color="black",ha="left", va="top", rotation=-20)
    n+=1
  axs[0].text(20,axs[0].get_ylim()[1],"TSS",color="r",ha="center",va="bottom")
  axs[0].text(100,axs[0].get_ylim()[1],"TES",color="b",ha="center",va="bottom")
  axs[1].set_ylabel("Methylation level")
  #fig.legend(loc=1, bbox_to_anchor=(1,1), bbox_transform=axs[0].transAxes,fancybox=True, framealpha=0.5)
  fig.savefig("methyl_function_region.pdf")

if __name__=="__main__":
  fire.Fire(plot_methyl)


