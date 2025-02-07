# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 16:38:48 2014
@author: Axel KournaK
This function convert matrices into triangular matrices 
to facilitate comparison with other signals.
"""
import os
from pylab import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import fire

def triangle(matrix,bdg,bed,start,end,res,triangled=True):
    start,end,res=map(int,[start,end,res]) # 起始终止分辨率int
    matrix_all=np.loadtxt(matrix) # 读入矩阵
    start_bin,end_bin=start/res+1,end/res+1 # 转到bin坐标，bin从1开始，原始位置从0开始，所以加1
    start_bin,end_bin=max(start_bin,1),min(end_bin,matrix_all.shape[1]) # 乱写的起始终止默认是全染色体
    bin_num=end_bin-start_bin+1 # bin总数

# 划分point，对应0.1/0.5/1/2/5/10M
    break_point=[100000,500000,1000000,2000000,5000000,10000000]
    nums=[abs((end-start)/i-5) for i in break_point] # 划分point数目和5的差值
    intersect=break_point[nums.index(min(nums))] # 取最小的point
    bin_point_num=intersect/res # 每个point的bin数目

    step_start=start_bin/bin_point_num+1 # point起始，从1开始
    step_end=end_bin/bin_point_num # point终止，最后的不要
    a,b=[],[] # x轴location,label
    if start % intersect == 0: # 整除
        a.append(0)
        b.append(str(float(start)/1000000)+"M") # label以M为单位
    bin_pos=step_start*bin_point_num-start_bin # 以bin为基本单位
    bin_pos_M=float(step_start*bin_point_num*res)/1000000
    for m in range(step_end-step_start+1): # 加1，range取到最后的前一个
        a.append(bin_point_num*m+bin_pos)
        pos=bin_pos_M+float(m*intersect/1000000)
        if pos%1==0: # 整数
            pos=int(pos)
        b.append(str(pos)+"M")

    gs=gridspec.GridSpec(3,1) # 区域
    matplotlib.rcParams.update({'font.size':8}) # 字体
    subplots_adjust(hspace=0,bottom=0.05,top=0.95) # 留白
    plt.figure(num=None,figsize=(6.52,5.20),dpi=150,facecolor="w",edgecolor="k") # 画布
    tri,vmean=matrix_all,np.log(np.nanmean(matrix_all)) # 矩阵，均值
    vmean_scale=float("1e"+"{:e}".format(vmean).split("e")[-1])
    if triangled: # 转三角形矩阵
        matrix_part=matrix_all[start_bin-1:,:end_bin]
        vmean=np.log(np.mean(matrix_part))
        tri=[]
        for i in range(0,matrix_part.shape[0],2):
            a=np.array([np.nan]*(i//2))
            a=np.concatenate((a,np.diagonal(matrix_part,offset=i),a))
            tri.append(a)
        tri=np.array(tri)
    np.savetxt("triangular.test.txt", tri, fmt='%.4f', delimiter=' ')

    ax1=plt.subplot(gs[:2,0]) # 区域
    im=ax1.imshow(np.log(tri[range(tri.shape[0]-1,-1,-1)]),cmap="OrRd",interpolation="none",vmin=vmean-0.25*vmean_scale,vmax=vmean+0.75*vmean_scale) # 热图
    plt.yticks([]) # 空y轴
    afx=plt.gca()
    afx.spines["top"].set_visible(False) # 抹除元素
    afx.spines["bottom"].set_visible(False)
    afx.spines["left"].set_visible(False)
    afx.spines["right"].set_visible(False)
    afx.tick_params(top=False,labelbottom=False,bottom=False,left=False,right=False)

    name=os.path.basename(matrix)
    abs_path=os.path.abspath(bdg)
    abs_odir=os.path.dirname(abs_path)

    start,end=map(str,[start,end])
    n=abs_odir+"/"+name+"-"+start+"-"+end+".triangular.insulation.png"
    n1=abs_odir+"/"+name+"-"+start+"-"+end+".triangular.insulation.pdf"

    ax2=plt.subplot(gs[2,0],sharex=ax1) # 区域
    plt.xticks(a,b) # x轴bin位置，实际位置label
    borders,borders2=[],[]
    with open(bdg) as bdg_handle:
        for bdg_line in bdg_handle:
            if bdg_line.startswith("track"):continue
            i=bdg_line.rstrip().split()
            bx=int(i[1])/res-start_bin # 以起点对应bin
            if i[3] == "NA":i[3]=0
            if bx<0 or bx>bin_num:continue
            by=float(i[3])
            borders.append(bx)
            borders2.append(by)
    ax2.set_xlim([0,bin_num]) # 范围
    ax2.plot(borders,borders2) # 线图
    
    with open(bed) as border:
        for line in border:
            if line.startswith("track"):continue
            i=line.rstrip().split()
            line=int(i[1])/res-start_bin
            if line<=0:continue
            plt.axvline(line,color="Gray",linewidth=.1,alpha=.5,linestyle="--") # 竖直线
    
    savefig(n,dpi=300,bbox_inches="tight")
    savefig(n1,format="pdf",dpi=300,bbox_inches="tight")
    close("all")

if __name__=="__main__":
    fire.Fire(triangle)


