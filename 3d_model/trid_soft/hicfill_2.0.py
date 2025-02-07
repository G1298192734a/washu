#!/public/frasergen/PUB/software/python/Python-3.7.10/bin/python3
from __future__ import division
import os
import sys
import scipy
from scipy import signal, stats
import numpy as np
import random
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from distfit import distfit


def find_flyout_multi_kernel(mat, 
                             marker=None,
                             out_f="find_flyout.png",
                             out_fig="find_flyout.heat.png"):
    find_flyout = find_mosaic
    plot_heat(mat, out_f=out_fig)
    diag = find_flyout(mat, size=3)
    diag_lst = []
    for size in range(3, 4, 1):
    #for size in [3, 9, 30]:
        kernel = np.concatenate([np.ones((size, size // 3)),
                                 -1 * np.ones((size, size // 3)),
                                 np.ones((size, size // 3))], axis=1)
        conv = signal.fftconvolve(mat, kernel, "same")
        diag_lst.append(np.diagonal(conv))
    diag = diag_lst[0]
    for diag_tmp in diag[1:]:
        diag = diag + diag_tmp
    plot_bar(diag, marker=marker, out_f=out_f)
    return diag

def find_flyout(mat, 
                marker=None, 
                size=3, 
                out_f="find_flyout.png", 
                out_fig="find_flyout.heat.png"):
    # marker [279, 341, 406, 439]
    plot_heat(mat, out_f=out_fig)
    kernel = np.concatenate([np.ones((size, size // 3)),
                             -1 * np.ones((size, size // 3)),
                             np.ones((size, size // 3))], axis=1)
    conv = signal.fftconvolve(mat, kernel, "same")
    diag = np.diagonal(conv)
    plot_bar(diag, marker=marker, out_f=out_f)
    return diag

def fit(array, distr="cauchy", out_prefix="dist"):
    dist = distfit(distr=distr)
    dist.fit_transform(array[array!=0])
    plt.clf()
    out_fig = out_prefix + ".best.png"
    dist.plot()
    plt.savefig(out_fig)    
    out_fig = out_prefix + ".flyout.png"
    plot_bar(array, thr=dist.model['CII_min_alpha'], out_f=out_fig)
    return dist

def find_mosaic(mat, 
                out_f="mosaic.png", 
                size=3, 
                marker=None, 
                out_fig="mosaic.heat.png"):
    mid = [random.choice([1, -1, 0]) for i in range(size * size//3)]    
    kernel = np.concatenate([np.ones((size, size // 3)),
                             np.reshape(mid, (size, size // 3)),
                             np.ones((size, size // 3))], axis=1)
    conv = signal.fftconvolve(mat, kernel, "same")
    diag = np.diagonal(conv)
    plot_bar(diag, marker=marker, out_f=out_f)    
    return diag

def plot_bar(array, out_f="bar_tmp.pdf", marker=None, thr=None):
    plt.clf()
    fig = plt.figure()
    axes = plt.gca()
    plt.plot(array)
    vmin = axes.get_ylim()[0]
    vmax = axes.get_ylim()[1]
    if marker:
        plt.vlines(marker, vmin, vmax, color="red")
    if thr:
        plt.hlines(thr, 0, len(array), color="orange")
    plt.savefig(out_f)

def plot_heat(mat, out_f="heat_tmp.pdf"):
    plt.clf()
    plt.imshow(mat, vmax=800)
    #plt.colorbar()
    plt.savefig(out_f)

def plot_together(mat, array, out_f="flyout.heat_and_bar.pdf", marker=None, thr=None):
    plt.clf()
    fig, axs = plt.subplots(2, 1, figsize=(15, 20), sharex='col', gridspec_kw={'height_ratios': [5, 1]})

    ax1 = axs[0]
    ax1.imshow(mat, vmax=80, cmap="OrRd")#, aspect='auto')#, interpolation='nearest')
    ax1.set_xlim(0, mat.shape[0])

    ax2 = axs[1]
    ax2.plot(array)
    vmin = ax2.get_ylim()[0]
    vmax = ax2.get_ylim()[1]
    if marker:
        ax2.vlines(marker, vmin, vmax, color="red")
    if thr:
        ax2.hlines(thr, 0, len(array), color="orange")
    ax2.set_xlim(0, mat.shape[0])
    ax2.set_ylim(vmin, vmax)
    
    plt.savefig(out_f)

def merge_region(ind):
    #print(ind)
    left = ind[0][0]
    right = ind[0][1]
    bound = []
    for i in ind[1:]:
        if i[0] - right <= 2:
            right = i[1]
        else:
            bound.append([left, right])
            left = i[0]
            right = i[1]
    bound.append([left, right])
    return bound

def choose_model(mat):
    score_lst = []
    for i in range(10):
        diag = find_mosaic(mat)
        dist = fit(diag)
        score_lst.append((dist, dist.model['score'], diag))
    return max(score_lst, key=lambda x: x[1])[0], max(score_lst, key=lambda x: x[1])[2]
    
def get_blank(mat, tolerate=3, ignore=10, bottom_params=0.3):
    #diag = find_flyout_multi_kernel(mat)
    #diag = find_mosaic(mat)
    #print(diag)
    #dist = fit(diag)
    dist, diag = choose_model(mat)
    all_ind = []
    #bottom = dist.model['CII_min_alpha'] + abs(dist.model['CII_min_alpha']) * 0.3
    bottom = dist.model['CII_min_alpha'] + abs(dist.model['CII_min_alpha']) * float(bottom_params)
    print("bottom:", bottom)
    for ind, val in enumerate(diag):
        if (val < bottom) and (ignore < ind < len(diag) - ignore):
            all_ind.append(ind)
    print("all_ind", all_ind)
    ind_lst = []
    try:
        for a, b in zip(all_ind, all_ind[1:]):
            if b - a <= tolerate:
                ind_lst.append([a, b])
        ind_lst = merge_region(ind_lst)
        done = [j for i in ind_lst for j in i]
    except:
        ind_lst = []
        done = []
    for i in all_ind:
        if i not in done:
            ind_lst.append(i)
    return ind_lst, list(diag), bottom

def fill(ind, mat, width=3):
    if type(ind) == list:
        l1 = ind[0]
        r1 = ind[1]
    else:
        l1 = ind
        r1 = ind
    l2 = l1 - width
    r2 = r1 + width
    #val = (mat[l2:l1] + mat[r1:r2]).mean(axis=0) / 2
    val = (np.mean(mat[l2:l1, :], axis=0) + np.mean(mat[r1:r2, :], axis=0)) / 2
    #val = np.nan_to_num(val)
    for i in range(l1, r1 + 1):
        mat[:, i] = val
        mat[i, :] = val
    return mat

def fill_blank(mat, ind_lst):
    for ind in ind_lst:
        mat = fill(ind, mat)
    return mat

def chrom_ind(chrom_len): # sorted by chrom ind
    bound = 0
    ind = [bound]
    for i in chrom_len:
        bound += i
        ind.append(bound)
    return ind

def split_heat(hic_heat, chrom_file, resolution=400000):    
    print("resolution:", resolution)
    chrom = pd.read_csv(chrom_file, header=None, sep="\t")
    chrom.columns = ["length"]
    chrom["length"] = chrom["length"] // resolution + 1 # resolution
    chrom_len = chrom["length"]
    
    left = chrom_ind(chrom_len)
    right = left[1:]
    heat_lst = []
    bound_lst = []
    for a, b in zip(left, right):
        heat_lst.append(hic_heat[a:b, a:b])
        bound_lst.append([a, b])
    return heat_lst, bound_lst  

def main(heat_file, chrom_file, bottom):
    heat = np.load(heat_file)
    heat_lst, bound_lst = split_heat(heat, chrom_file, resolution=RESOLUTION)
    heat_final = heat.copy()
    count = 0
    chrom_size = open(chrom_file).readlines()
    chrom_size = [i.strip() for i in chrom_size]
    f1 = open(chrom_size_file, "w+")
    f2 = open(chrom_name_file, "w+")
    for mat, bound, size in zip(heat_lst, bound_lst, chrom_size):
        count += 1
        a, b  = bound
        ind_lst, score, bottom = get_blank(mat, ignore=1, bottom_params=bottom) #IGNORE = 0
        marker = []
        for i in ind_lst:
            if type(i) == list:
                for j in i:
                    marker.append(j)
            else:
                marker.append(i)
        #plot_together(mat, score, out_f="heat_with_bar.before.{}.png".format(count), marker=marker, thr=bottom)
        #plot_together(mat, score, out_f="heat_with_bar.before.{}.pdf".format(count), marker=marker, thr=bottom)
        mat = fill_blank(mat, ind_lst)
        heat_final[a:b, a:b] = mat
        #plot_together(mat, score, out_f="heat_with_bar.after.{}.png".format(count), marker=marker, thr=bottom)
        #plot_together(mat, score, out_f="heat_with_bar.after.{}.pdf".format(count), marker=marker, thr=bottom)
        f1.write(size + "\n")
        f2.write(str(count) + "\n")
    f1.close()
    f2.close()
    return heat_final

if __name__ == "__main__":
    heat_file       = sys.argv[1]
    chrom_file      = sys.argv[2]
    RESOLUTION      = int(sys.argv[3])
    out_f           = sys.argv[4]
    chrom_size_file = sys.argv[5]
    chrom_name_file = sys.argv[6]
    bottom          = sys.argv[7]
    heat            = main(heat_file, chrom_file, bottom)
    np.save(out_f, heat)        

