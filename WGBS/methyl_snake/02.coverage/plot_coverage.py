import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import fire

def read_csv_chunk(filePath):
  dfs=pd.read_table(filePath,chunksize=1000000,dtype="int",usecols=[2],header=None,engine="c",na_filter=False,low_memory=False)
  depth_all,genome_len=[], 0
  for i in dfs:
    depth=i[2].value_counts()
    genome_len+=depth.sum()
    depth_all.append(depth[depth.index<=50])
  return genome_len,pd.concat(depth_all,axis=1,join="outer").sum(axis=1)

def plot_coverage(depth):
  genome_len,plot_data=read_csv_chunk(depth)
  plot_data=plot_data.sort_index()

  fig, axs = plt.subplots(1, 1)
  spl=make_interp_spline(plot_data.index,plot_data.tolist()/genome_len)
  x_smooth=np.linspace(0,50,300)
  y_smooth=spl(x_smooth)
  axs.plot(x_smooth,y_smooth,label="distribution")
  axs.set_ylabel("distribution")
  axs.set_title("Coverage")

  def accumu():
    a,b=0,genome_len
    while True:
      yield b
      a,b=a+1,b-plot_data[a]
  acm=accumu()
  spl1=make_interp_spline(plot_data.index,[next(acm)/genome_len for i in range(51)])
  y_smooth1=spl1(x_smooth)
  axs1=axs.twinx()
  axs1.plot(x_smooth,y_smooth1,c="red",label="accumulation")
  axs1.set_ylabel("accumulation")
  fig.legend(loc=1, bbox_to_anchor=(1,1), bbox_transform=axs.transAxes,fancybox=True, framealpha=0.5)
  plt.savefig("coverage_distribution.pdf")

if __name__=="__main__":
  fire.Fire(plot_coverage)

