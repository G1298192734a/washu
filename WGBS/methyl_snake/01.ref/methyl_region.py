import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig as pb
from pybedtools import BedTool as bt
from scipy.interpolate import make_interp_spline
from itertools import chain
import fire

def get_methyl(bw,*beds):
  bw=pb.open(bw)
  methyl_avg=[*chain(*(pd.DataFrame([
    [bw.stats(j.chrom,p,q,exact=True)[0] for p,q in zip(np.linspace(j.start,j.end,21,dtype=int),np.linspace(j.start,j.end,21,dtype=int)[1:])]
    for j in i if j.length>21
  ]).mean().tolist() for i in map(bt,beds)))]
  x=np.arange(20*len(beds))
  spl=make_interp_spline(x,methyl_avg)
  x_smooth=np.linspace(x.min(),x.max(),300)
  y_smooth=spl(x_smooth)
  plt.plot(x,methyl_avg,"o",label="Original Data")
  plt.plot(x_smooth,y_smooth,label="Smoothed Line")
  plt.legend()
  plt.savefig("test.pdf")

if __name__=="__main__":
  fire.Fire(get_methyl)


