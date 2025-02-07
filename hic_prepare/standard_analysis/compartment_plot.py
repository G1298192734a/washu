import argparse, os,textwrap
import cooler 
import numpy as np
import pandas as pd
from cooltools.lib import numutils
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
matplotlib.rc("font", size=24)
import warnings
warnings.filterwarnings("ignore")

__author__  = "lurui@frasergen.com"
__date__    = "2023.03.16"
__version__ = "v0.1"

# date: 2024.01.17
# change: use argparse instead of fire, add raw matrix support 

def get_raw_matrix(incool, chrom):
	mycool = cooler.Cooler(incool)
	cool_matrix = mycool.matrix(balance=True).fetch(chrom)
	cool_matrix = np.array(cool_matrix)
	cool_matrix[~np.isfinite(cool_matrix)] = 0
	return cool_matrix

def get_oe_matrix(incool, chrom):
	""" [NOTE] this observed over expected calculation come from:
	           cooltools.api.eigdecomp.cis_eig 
	"""
	mycool = cooler.Cooler(incool)
	cool_matrix = mycool.matrix(balance=True).fetch(chrom)
	cool_matrix = np.array(cool_matrix)
	cool_matrix[~np.isfinite(cool_matrix)] = 0
	mask = cool_matrix.sum(axis=0) > 0
	return numutils.observed_over_expected(cool_matrix, mask)[0]

def one_sample(matrix, resolution, chrom, outdir, name1, pc1_table1, min_value, max_value):
	max_len = matrix.shape[0]
	fruitpunch = sns.blend_palette(["white", "red"], as_cmap=True)

	fig=plt.figure(figsize=(12, 14))
	matplotlib.rcParams["font.size"] = 18
	gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2,15], hspace=0)

	plt.subplot(gs[0, 0])
	frame = plt.gca()
	frame.spines['top'].set_visible(False)
	frame.spines['right'].set_visible(False)
	frame.spines['bottom'].set_position(('data',0))
	frame.axes.get_xaxis().set_visible(False)
	frame.set_xlim(0,max_len)
	frame.set_ylim(-2, 2)
	frame.set_yticks([-2, 0, 2])
	frame.set_yticklabels(["", "0", "2"])
	frame.bar(data=pc1_table1, x="locus", height="PC1", color="color", width=1)

	plt.subplot(gs[1, 0])
	ax=plt.gca()
	image=ax.imshow(matrix, interpolation='nearest', cmap=fruitpunch, rasterized=True, vmin=min_value, vmax=max_value)
	ax.set_xlim(0, max_len)
	ax.set_ylim(max_len, 0)
	ax.plot([0, max_len], [0, max_len], color="black", ls="--", alpha=0.8, lw=1)
	ticks_lab = (ax.get_xticks() * resolution).astype(int)
	if np.all(ticks_lab % 1000000 == 0):
		ax.set_xticklabels([f"{lab // 1000000}Mb" for lab in ticks_lab])
		ax.set_yticklabels([f"{lab // 1000000}Mb" for lab in ticks_lab])

	cax = fig.add_axes([ax.get_position().x1+0.05,ax.get_position().y0,0.03,ax.get_position().height])
	plt.colorbar(image, cax=cax, shrink=0.9, fraction=0.6, aspect=30)

	plt.savefig(os.path.join(outdir, f"{name1}.{chrom}.png"), dpi=300, bbox_inches="tight")
	plt.savefig(os.path.join(outdir, f"{name1}.{chrom}.pdf"), bbox_inches="tight",)
	
def two_sample(matrix, resolution, chrom, outdir, name1, pc1_table1, min_value, max_value, name2, pc1_table2):
	max_len = matrix.shape[0]
	fruitpunch = sns.blend_palette(["white", "red"], as_cmap=True)

	fig=plt.figure(figsize=(12, 14))
	matplotlib.rcParams["font.size"] = 18
	#gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[2, 2, 20], width_ratios=[15, 1], wspace=0, hspace=0)
	gs = gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[2, 2, 15], width_ratios=[9, 1], wspace=0, hspace=0.02)

	control_cptmt = plt.subplot(gs[0, 0])
	frame = plt.gca()
	frame.spines['top'].set_visible(False)
	frame.spines['right'].set_visible(False)
	frame.spines['bottom'].set_position(('data',0))
	frame.axes.get_xaxis().set_visible(False)
	frame.set_xlim(0,max_len)
	frame.set_ylim(-2, 2)
	frame.bar(data=pc1_table1, x="locus", height="PC1", color="color", width=1)
	frame.set_yticks([-2, 0, 2])
	frame.set_yticklabels(["", "0", "2"])

	axes = plt.subplot(gs[0, 1])
	axes.text(0.1, 0.5, textwrap.fill(name1,9), fontsize=24, ha="left", va="center", transform=plt.gca().transAxes)	
	axes.set_xticks([])
	axes.set_yticks([])
	for place in ["top", "right", "bottom", "left"]:
		axes.spines[place].set_visible(False)

	case_cptmt = plt.subplot(gs[1, 0]) 
	frame = plt.gca()
	frame.spines['top'].set_visible(False)
	frame.spines['right'].set_visible(False)
	frame.spines['bottom'].set_position(('data',0))
	frame.axes.get_xaxis().set_visible(False)
	frame.set_xlim(0,max_len)
	frame.set_ylim(-2, 2)
	frame.bar(data=pc1_table2, x="locus", height="PC1", color="color", width=1)
	frame.set_yticks([-2, 0, 2])
	frame.set_yticklabels(["", "0", "2"])

	axes = plt.subplot(gs[1, 1])
	axes.text(0.1, 0.5, textwrap.fill(name2,9), fontsize=24, ha="left", va="center", transform=plt.gca().transAxes)		
	axes.set_xticks([])
	axes.set_yticks([])
	for place in ["top", "right", "bottom", "left"]:
		axes.spines[place].set_visible(False)

	plt.subplot(gs[2, 0])
	axes = plt.gca()
	image = axes.imshow(matrix, interpolation='nearest', cmap=fruitpunch, rasterized=True, vmin=min_value, vmax=max_value)
	axes.text(0.6, 0.8, textwrap.fill(name1,9), fontsize=30, transform=plt.gca().transAxes)
	axes.text(0.2, 0.2, textwrap.fill(name2,9), fontsize=30, transform=plt.gca().transAxes)
	axes.set_xlim(0, max_len)
	axes.set_ylim(max_len, 0)
	axes.plot([0, max_len], [0, max_len], color="black", ls="--", alpha=0.8, lw=1)
	ticks_lab = (axes.get_xticks() * resolution).astype(int)
	if np.all(ticks_lab % 1000000 == 0):
	    axes.set_xticklabels([f"{lab // 1000000}Mb" for lab in ticks_lab])
	    axes.set_yticklabels([f"{lab // 1000000}Mb" for lab in ticks_lab])
	
	#divider = make_axes_locatable(axes)
	#cax = divider.append_axes("right", size="3%", pad=0.3)
	cax = fig.add_axes([axes.get_position().x1+0.05,axes.get_position().y0,0.03,axes.get_position().height])
	plt.colorbar(image, cax=cax, shrink=0.9, fraction=0.6, aspect=30)
	
	plt.savefig(os.path.join(outdir, f"{name2}.minus.{name1}.{chrom}.png"), dpi=300, bbox_inches="tight")
	plt.savefig(os.path.join(outdir, f"{name2}.minus.{name1}.{chrom}.pdf"), bbox_inches="tight")

def generate_heatmap(matrix, resolution, chrom, outdir, name1, pc1_table1, 
					min_value, max_value, name2=None, pc1_table2=None):
	if name2:
		two_sample(matrix, resolution, chrom, outdir, name1, pc1_table1, min_value, max_value, name2, pc1_table2)
	else: 
		one_sample(matrix, resolution, chrom, outdir, name1, pc1_table1, min_value, max_value)
		
def get_cptmt_pc1(cptmt_gcgene, chrom, pc1_column, resolution):
	cptmt_pc1 = pd.read_csv(cptmt_gcgene, usecols=(0, 1, 2, pc1_column-1), sep="\t",
					names=("chrom", "start", "end", "PC1"))
	cptmt_pc1 = cptmt_pc1[cptmt_pc1.chrom==chrom]
	cptmt_pc1["locus"] = (cptmt_pc1.start + cptmt_pc1.end) // 2
	cptmt_pc1["color"] = 'steelblue'
	cptmt_pc1.loc[cptmt_pc1.PC1==0, "color"] = "grey"
	cptmt_pc1.loc[cptmt_pc1.PC1<0, "color"] = "firebrick"
	cptmt_pc1.locus = cptmt_pc1.locus / resolution
	return cptmt_pc1[["locus", "PC1", "color"]]
	
def main(myargs):
	resolution = cooler.Cooler(myargs.control_cool).binsize
	
	if not os.path.exists(myargs.outdir): os.makedirs(myargs.outdir)
	
	for line in open(myargs.chromsize):
		chrom = line.strip().split()[0]
		if myargs.treat_name and myargs.treat_cool: 
			if myargs.matrix_type == "OE":
				control_df = get_oe_matrix(myargs.control_cool, chrom)
				treat_df = get_oe_matrix(myargs.treat_cool, chrom)
			elif myargs.matrix_type == "raw":
				control_df = get_raw_matrix(myargs.control_cool, chrom)
				treat_df = get_raw_matrix(myargs.treat_cool, chrom)
			control_pc1 = get_cptmt_pc1(myargs.control_cptmt_gcgene, chrom, myargs.pc1_column, resolution)
			treat_pc1 = get_cptmt_pc1(myargs.treat_cptmt_gcgene, chrom, myargs.pc1_column, resolution)
			#generate_heatmap(treat_oe, resolution, outdir, treat_name, treat_pc1)
		
			control_up = np.triu(np.array(control_df), k=1)
			treat_low = np.tril(np.array(treat_df), k=-1)
			combine = control_up + treat_low
			generate_heatmap(combine, resolution, chrom, myargs.outdir, myargs.control_name, 
						control_pc1, myargs.min_value, myargs.max_value, myargs.treat_name, treat_pc1)	
		else:
			if myargs.matrix_type == "OE":
				control_df = get_oe_matrix(myargs.control_cool, chrom)
			elif myargs.matrix_type == "raw":
				control_df = get_raw_matrix(myargs.control_cool, chrom)
			control_pc1 = get_cptmt_pc1(myargs.control_cptmt_gcgene, chrom, myargs.pc1_column, resolution)
			generate_heatmap(control_df, resolution, chrom, myargs.outdir, myargs.control_name, control_pc1, 
						myargs.min_value, myargs.max_value)

def get_myargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--chromsize", dest="chromsize", required=True, 
			help="reference chrom size file [required]")
	parser.add_argument("-n1", "--control_name", dest="control_name", required=True, 
			help="control sample name [required]")
	parser.add_argument("-c1", "--control_cool", dest="control_cool", required=True, 
			help="control sample cool file, eg. WT.100000.cool [required]")
	parser.add_argument("-g1", "--control_cptmt_gcgene", dest="control_cptmt_gcgene", required=True, 
			help="control sample compartment GC gene table, eg.WT_100000.gc_gene.bed [requried]")
	parser.add_argument("-n2", "--treat_name", dest="treat_name", default=None, 
			help="treat sample name [default:None]")
	parser.add_argument("-c2", "--treat_cool", dest="treat_cool", default=None, 
			help="treat sample cool file, eg. MUT.100000.cool [default:None]")
	parser.add_argument("-g2", "--treat_cptmt_gcgene", dest="treat_cptmt_gcgene", default=None, 
			help="treat sample compartment GC gene table, eg. MUT_100000.gc_gene.bed [default:None]")
	parser.add_argument("-m", "--matrix_type", dest="matrix_type", default="OE", choices=["OE", "raw"], 
			help="which heatmap to draw, choose from OE matrix or raw matrix [default:OE]")
	parser.add_argument("-x", "--min_value", dest="min_value", default=0, type=float, 
			help="min value for HiC interaction heatmap [default:0]")
	parser.add_argument("-y", "--max_value", dest="max_value", default=2, type=float, 
			help="max value for HiC interaction heatmap [default:2]")
	parser.add_argument("-c", "--pc1_column", dest="pc1_column", default=4, type=int, 
			help="which column record PC1 value in GC gene table [default:4]")
	parser.add_argument("-o", "--outdir", dest="outdir", default="./", 
			help="output dir [default:./]")
	return parser.parse_args()

if __name__ == "__main__":
	main(get_myargs())
