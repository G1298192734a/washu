import argparse, gzip
import cooler, os
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
matplotlib.rc("font", size=14, family="serif")

def parse_genome_cool(incool, outmatrix):
	mycool = cooler.Cooler(incool)
	bin_table = mycool.bins()[:]
	bin_table["bin_count"] = (bin_table.index + 1).astype(str)
	bin_table["location"] = bin_table.chrom.astype(str) + ":" + bin_table.start.astype(str) +"-"+ bin_table.end.astype(str)
	bin_table["header"] = "HiC_" + bin_table.bin_count + "|XXX|" + bin_table.location

	interaction = mycool.matrix(as_pixels=False, balance=True)
	with gzip.open(outmatrix, "wt") as outbuff:
		print("\t" + "\t".join(bin_table.header), file=outbuff)

		for i, line in enumerate(interaction):
			line = np.array(line)
			line[np.isnan(line)] = 0
			print(bin_table.header.values[i], "\t".join(map(str, line[0])), sep="\t", file=outbuff)
	print("[NOTE] parse cool file: done !")
	cmd = f"/public/frasergen/PUB/software/perl/perl-5.34.0/bin/perl -I /work/frasergen/3D/pipeline/Interactome/hic3_workflow/src/cworld-dekker/cworld-dekker/lib /work/frasergen/3D/pipeline/Interactome/hic3_workflow/src/cworld-dekker/cworld-dekker/scripts/perl/matrix2loess.pl -i {outmatrix} -o {outmatrix.replace('.matrix.gz','')} --ca 0.005 --ez && echo [NOTE] generate zscore matrix: done !"
	print(cmd)
	os.system(cmd)

def parse_chrom_cool(incool, outmatrix):
	mycool = cooler.Cooler(incool)
	bin_table = mycool.bins()[:]
	bin_table["bin_count"] = (bin_table.index + 1).astype(str)
	bin_table["location"] = bin_table.chrom.astype(str) + ":" + bin_table.start.astype(str) +"-"+ bin_table.end.astype(str)
	bin_table["header"] = "HiC_" + bin_table.bin_count + "|XXX|" + bin_table.location

	chrom_list = bin_table.chrom.unique()
	interaction = mycool.matrix(as_pixels=False, balance=True)
	for chrom in chrom_list:
		chrom_bin = bin_table[bin_table.chrom==chrom]
		chrom_matrix = pd.DataFrame(interaction.fetch(chrom), index=list(chrom_bin.header), 
						columns=list(chrom_bin.header)).fillna(0)
		chrom_outfile = str(outmatrix).replace(".matrix.gz", f".{chrom}.matrix.gz")
		chrom_matrix.to_csv(chrom_outfile, sep="\t", compression="gzip")
		print(f"[NOTE] parse {chrom} cool file: done !")
		cmd = f"/public/frasergen/PUB/software/perl/perl-5.34.0/bin/perl -I /work/frasergen/3D/pipeline/Interactome/hic3_workflow/src/cworld-dekker/cworld-dekker/lib /work/frasergen/3D/pipeline/Interactome/hic3_workflow/src/cworld-dekker/cworld-dekker/scripts/perl/matrix2loess.pl -i {chrom_outfile} -o {chrom_outfile.replace('.matrix.gz','')} --ca 0.005 --ez && echo [NOTE] generate {chrom} zscore matrix: done !"
		os.system(cmd)

def plot_matrix(matrix_type, delta_matrix, binsize, outpfix, chrom_len=None):
	plt.figure(figsize=(10, 10))
	heatmap = plt.imshow(delta_matrix, interpolation='nearest', cmap='seismic',vmin=-10, vmax=10)
	plt.colorbar(heatmap, shrink=0.8)
	
	if matrix_type == "genome":
		chrom_bin = np.ceil(chrom_len/binsize).astype(int)
		max_bin = chrom_bin.sum()
		chrom_cums = chrom_bin.cumsum()[:-1]
		plt.vlines(x=chrom_cums, ymin=0, ymax=max_bin, color="grey", 
					linestyles="dashed", linewidth=0.6, alpha=0.75)
		plt.hlines(y=chrom_cums, xmin=0, xmax=max_bin, color="grey", 
					linestyles="dashed", linewidth=0.6, alpha=0.75)

	plt.xlim(0, delta_matrix.shape[0])
	plt.ylim(delta_matrix.shape[1], 0)
	axes = plt.gca()

	old_xticks = axes.get_xticks()[:-1]
	new_xticks = [int(i) if i>1 else i for i in old_xticks * binsize / 1000000]
	plt.xticks(ticks=old_xticks, labels=new_xticks)
	
	old_yticks = axes.get_yticks()[:-1]
	new_yticks = [int(i) if i>1 else i for i in old_yticks * binsize / 1000000]
	plt.yticks(ticks=old_yticks, labels=new_yticks)
	plt.xlabel("Genomic position (Mb)")
	plt.ylabel("Genomic position (Mb)")
	plt.savefig(f"{outpfix}.png", dpi=300, bbox_inches="tight")
	plt.savefig(f"{outpfix}.pdf", bbox_inches="tight")
	plt.close()

def genome_diff_matrix(myargs, outdir, chrom_len):
	wt_zscore = pd.read_csv(myargs.wt_mat, sep="\t", index_col=0, header=0, comment="#")
	case_zscore_fn = outdir.joinpath(f"{myargs.case_name}.zScore.matrix.gz")
	case_zscore = pd.read_csv(myargs.case_mat, sep="\t", index_col=0, header=0, comment="#")
	delta_matrix = case_zscore - wt_zscore
	out_matrix = outdir.joinpath(f"{myargs.case_name}.minus.{myargs.wt_name}.matrix.gz")
	delta_matrix.to_csv(out_matrix, sep="\t", index=True, header=True)
	print("[NOTE] generate differential matrix: done !")

	fig_prefix = outdir.joinpath(f"{myargs.case_name}.minus.{myargs.wt_name}.genome")
	plot_matrix(myargs.matrix, delta_matrix, int(myargs.res), fig_prefix, chrom_len)
	print("[NOTE] generate differential heatmap: done !")

def chrom_diff_matrix(myargs, outdir, chrom_list):
	for chrom in chrom_list:
		wt_zscore = pd.read_csv(Path(myargs.wt_mat).joinpath(f"{myargs.wt_name}.{chrom}.zScore.matrix.gz"), sep="\t", index_col=0, header=0, comment="#")
		case_zscore = pd.read_csv(Path(myargs.case_mat).joinpath(f"{myargs.case_name}.{chrom}.zScore.matrix.gz"), sep="\t", index_col=0, header=0, comment="#")
		delta_matrix = case_zscore - wt_zscore
		out_matrix = outdir.joinpath(f"{myargs.case_name}.minus.{myargs.wt_name}.{chrom}.matrix.gz")
		delta_matrix.to_csv(out_matrix, sep="\t", index=True, header=True)
		print(f"[NOTE] generate {chrom} differential matrix: done !")
	
		fig_prefix = outdir.joinpath(f"{myargs.case_name}.minus.{myargs.wt_name}.{chrom}")
		plot_matrix(myargs.matrix, delta_matrix, int(myargs.res), fig_prefix)
		print(f"[NOTE] generate {chrom} differential heatmap: done !")

def get_myargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--matrix", choices={"genome", "chrom"}, required=True, 
				help="genome or chrom differential matrix to generate, choose from `genome' or `chrom' [required]")
	parser.add_argument("-w", "--wt_mat", dest="wt_mat", required=True, 
			help="wildtype or control sample zscore matrix file [required]")
	parser.add_argument("-c", "--case_mat", dest="case_mat", required=True, 
			help="case or treat sample zscore matrix file [required]")
	parser.add_argument("-wn", "--wt_name", dest="wt_name", required=True, 
				help="wildtype sample name [required]")
	parser.add_argument("-cn", "--case_name", dest="case_name", required=True, 
				help="case sample name [required]")
	parser.add_argument("-o", "--outdir", dest="outdir", default="./", 
				help="output dir [defult:./]")
	parser.add_argument("-g", "--csize", dest="csize", required=True, 
				help="genome.chromsize [required]")
	parser.add_argument("-r", "--res", dest="res", required=True, 
				help="resolution [required]")
	return parser.parse_args()

if __name__ == "__main__":
	myargs = get_myargs()
	outdir = Path(myargs.outdir).joinpath(f"{myargs.case_name}.minus.{myargs.wt_name}")
	outdir.mkdir(parents=True, exist_ok=True)

	csize = pd.read_table(myargs.csize,header=None,engine="c")
	chrom_list,chrom_len = csize[0],csize[1]
	
	if myargs.matrix == "genome":
		genome_diff_matrix(myargs, outdir, chrom_len)
	elif myargs.matrix == "chrom":
		chrom_diff_matrix(myargs, outdir, chrom_list)
