# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:39:30 2015
@author: Axel KournaK
2 DI + SCALOGRAM
"""
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from pylab import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import fire

def ice(A,threshold):
    n,a=A.shape[0],A.sum(axis=0)
    for i in range(51):
        b=np.sum(A,axis=0)
        A/=b*b.reshape(n,1)
        A[np.isnan(A)]=0
        A[a<threshold]=0
        A[:,a<threshold]=0
    return A

def scn(A,threshold):
    a=A.sum(axis=0)
    A/=a
    A/=A.sum(axis=0)
    A[np.isnan(A)]=0
    A[a<threshold]=0
    A[:,a<threshold]=0
    return A

def comp(A,nw,circ):
    n=A.shape[0]
    a=A if circ else np.zeros_like(A)
    A=np.concatenate((a,A,a),axis=1)
    B=np.array([np.diagonal(A,offset=i) for i in range(n-nw+1,n+nw)]).T
    C=np.array([B[:,nw-i-1:nw+i].sum(axis=1) for i in range(nw)])
    return C

def __norm__(matrix,range_bin):
    matscn=ice(matrix,100)
    matscn=scn(matscn,0)
    comp_scales1=comp(matscn,range_bin)
    np.savetxt("compact.txt", comp_scales1, fmt="%.4f",delimiter=" ")
    return comp_scales1

def compact(matrix,res,output,rbin=400,norm=True,circ=True):
    res,range_bin=int(res)/1000,int(rbin)
    data=np.loadtxt(matrix)
    comp_scales1=__norm__(data,range_bin,circ) if norm else data
    gs=gridspec.GridSpec(1,1)
    matplotlib.rcParams.update({'font.size':8})
    subplots_adjust(hspace=0)
    plt.figure(num=None,figsize=(10,3),dpi=800,facecolor="w",edgecolor="k")
    ax3=plt.subplot(gs[0])
    xlim,ylim=data.shape[1],range_bin
    ax3.contourf(comp_scales1,vmin=0.0,vmax=1.0,levels=[0,.15,.30,.45,.6,.75,1.0],extent=[0,xlim,0,ylim],cmap="rainbow")
    im=contourf(comp_scales1,vmin=0.0,vmax=1.0,levels=[0,.15,.30,.45,.6,.75,1.0],extent=[0,xlim,0,ylim],cmap="rainbow")
    ax3.set_ylabel("Scales (in kb)")
    tick_locs=range(0,ylim,ylim/10)
    tick_lbls=(np.array(tick_locs)*res).tolist()
    plt.yticks(tick_locs,tick_lbls,fontsize=10)
    xlim_len=len(str(xlim*res))-1
    tick_locs=range(0,xlim,10**xlim_len/res)
    tick_lbls=(np.array(tick_locs)*res).tolist()
    plt.xticks(tick_locs,tick_lbls,fontsize=10)
    bounds=[0.,1.0]
    cb1=colorbar(im,shrink=.2,orientation="horizontal",ticks=bounds,spacing="proportional")
    plt.savefig(output+"_DOM"+".svg",dpi=600,format="svg")
    plt.savefig(output+"_DOM"+".eps",dpi=600,format="eps")
    plt.savefig(output+"_DOM"+".jpeg",dpi=600,format="jpeg")
    close("all")

if __name__=="__main__":
    fire.Fire(compact)


