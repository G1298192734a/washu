import argparse
import cooler
from pathlib import Path
import numpy as np 
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
matplotlib.rc("font", size=14, family="serif")
import seaborn as sns


def plot_matrix(matrix_type, in_matrix, binsize, outpfix, vmin, vmax, 
				base_color, chrom_len=None):
	plt.figure(figsize=(10, 10))
	fruitpunch = sns.blend_palette(["white", base_color], as_cmap=True)

	if vmin >= 0:
		heatmap = plt.imshow(in_matrix, interpolation='nearest', 
						cmap=fruitpunch, vmin=vmin, vmax=vmax)
	else: 
		heatmap = plt.imshow(in_matrix, interpolation='nearest',
						cmap="bwr", vmin=vmin, vmax=vmax)
	plt.colorbar(heatmap, shrink=0.8)
	
	if matrix_type == "genome":
		chrom_bin = np.ceil(chrom_len/binsize).astype(int)
		max_bin = chrom_bin.sum()
		chrom_cums = chrom_bin.cumsum()[:-1]
		plt.vlines(x=chrom_cums, ymin=0, ymax=max_bin, color="grey", 
					linestyles="dashed", linewidth=0.6, alpha=0.75)
		plt.hlines(y=chrom_cums, xmin=0, xmax=max_bin, color="grey", 
					linestyles="dashed", linewidth=0.6, alpha=0.75)

	plt.xlim(0, in_matrix.shape[0])
	plt.ylim(in_matrix.shape[1], 0)
	axes = plt.gca()

	old_xticks = axes.get_xticks()[:-1]
	new_xticks = [ int(i) if i>1 else i for i in old_xticks * binsize / 1000000]
	plt.xticks(ticks=old_xticks, labels=new_xticks)
	
	old_yticks = axes.get_yticks()[:-1]
	new_yticks = [ int(i) if i>1 else i for i in old_yticks * binsize / 1000000]
	plt.yticks(ticks=old_yticks, labels=new_yticks)
	plt.xlabel("Genomic position (Mb)")
	plt.ylabel("Genomic position (Mb)")
	plt.savefig(f"{outpfix}.png", dpi=300, bbox_inches="tight")
	plt.savefig(f"{outpfix}.pdf", bbox_inches="tight")
	plt.close()

def main(myargs):
	Path(myargs.outpfix).parent.mkdir(parents=True, exist_ok=True)

	mycool = cooler.Cooler(myargs.incool)
	chrom_len = mycool.chromsizes
	binsize = mycool.binsize
	if myargs.region == "genome":
		if myargs.rawmatrix:
			genome_cool = mycool.matrix(balance=False)[:]
		else:
			genome_cool = mycool.matrix(balance=True)[:]
		plot_matrix("genome", genome_cool, binsize, f"{myargs.outpfix}.genome", 
				myargs.vmin, myargs.vmax, myargs.basecolor, chrom_len)
	else:
		for chrom in mycool.chromnames:
			if myargs.rawmatrix:
				chrom_cool = mycool.matrix(balance=False).fetch(chrom)
			else:
				chrom_cool = mycool.matrix(balance=True).fetch(chrom)
			plot_matrix("chrom", chrom_cool, binsize, f"{myargs.outpfix}.{chrom}",
				myargs.vmin, myargs.vmax, myargs.basecolor, chrom_len)

def get_myargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--incool", required=True, 
			help="input Hi-C matrix cool file [required]")
	parser.add_argument("-g", "--region", choices={"genome", "chrom"}, 
			help="generate genome heatmap or chromosome heatmap [required]")
	parser.add_argument("-r", "--rawmatrix", action="store_true", default=False, 
			help="use cool raw matrix or balanced matrix(default)")
	parser.add_argument("--vmin", default=0, type=float, 
			help="min value for heatmap colorbar [default:0]")
	parser.add_argument("--vmax", default=0.03, type=float, 
			help="max value for heatmap colorbar [default:0.03]")
	parser.add_argument("--basecolor", default="red", 
			help="base color for heatmap [default: red]")
	parser.add_argument("-o", "--outpfix", default="sampleid", 
			help="output prefix [default:sampleid]")
	return parser.parse_args()
	

if __name__ == "__main__":
	main(get_myargs())

