import argparse, os
import bioframe
import pandas as pd
import numpy as np
from collections import OrderedDict
from itertools import permutations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns 
plt.rc("font", family='serif', size=10)

# 20240117. change colormap to Reds

def get_reproduce(dataset_pairs, dataset_contact, genome_bin, genomedisco, outdir):
	for biorep1, biorep2 in dataset_pairs:
		vs_name = f"{biorep1}.vs.{biorep2}"
		out_vsdir = os.path.join(outdir, vs_name)
		try: os.stat(out_vsdir)
		except: os.makedirs(out_vsdir)

		contact_file = os.path.join(out_vsdir, f"{vs_name}.contact")
		with open(contact_file, "w") as contact_buff:
			print(biorep1, dataset_contact[biorep1], sep="\t", file=contact_buff)	
			print(biorep2, dataset_contact[biorep2], sep="\t", file=contact_buff)

		pair_file = os.path.join(out_vsdir, f"{vs_name}.pair")
		with open(pair_file, "w") as pair_buff:
			print(biorep1, biorep2, sep="\t", file=pair_buff)
	
		shell = os.path.join(out_vsdir, f"{vs_name}.sh")	
		with open(shell, "w") as outbuff:
			cmd = (f"{genomedisco} run_all --metadata_samples {contact_file} "
				f"--metadata_pairs {pair_file} --bins {genome_bin} --outdir {out_vsdir}")
			print("#!/bin/bash", file=outbuff)
			print(cmd, file=outbuff)
			os.system(cmd)

def get_plot(dataset_pairs, dataset_list, outdir,labels,label_lens):
	total_cor = dict()
	xcax=len(max(dataset_list,key=len))*0.013 if len(dataset_list)>10 else 0.03
	for biorep1, biorep2 in dataset_pairs:
		vs_name = f"{biorep1}.vs.{biorep2}"
		corr_file = os.path.join(outdir, vs_name, "results", "reproducibility", "GenomeDISCO", f"{vs_name}.txt")
		corr_df = pd.read_csv(corr_file, sep="\t", names=("biorep1", "biorep2", "chrom", "correlation"))	
		corr_value = corr_df.correlation.mean()
		total_cor[(biorep1, biorep2)] = corr_value
		total_cor[(biorep2, biorep1)] = corr_value

	corr_list = list()
	for i in dataset_list: 
		temp_list = list()
		for j in dataset_list:
			if i == j: 
				temp_list.append(1)
			else: 
				temp_list.append(total_cor[(i, j)])
		corr_list.append(temp_list)

	mydf = pd.DataFrame(corr_list, index=dataset_list, columns=dataset_list)
	fig=plt.figure(figsize=(12, 10))
	#fruitpunch = sns.blend_palette(["white", "red"], as_cmap=True)
	#sns.heatmap(mydf, cmap=fruitpunch, annot=True, fmt=".3f", annot_kws={"fontsize": 18},
	#sns.clustermap(mydf, cmap="Reds", annot=True, fmt=".3f", annot_kws={"fontsize": 25/np.sqrt(mydf.shape[0])},
	ax=sns.heatmap(mydf, cmap="Reds", annot=True, fmt=".3f", annot_kws={"fontsize": 25/np.sqrt(mydf.shape[0])},
			xticklabels=True, yticklabels=True)
	pos = 0
	cax=fig.add_axes([ax.get_position().x0-xcax,ax.get_position().y0,0.005,ax.get_position().height])
	cax.set_ylim(0,sum(label_lens))
	cax.set_xticks([])
	cax.set_yticks([])
	for place in ["top", "right", "bottom", "left"]:
		cax.spines[place].set_visible(False)
	for i in label_lens[::-1]:
		cax.vlines(0,pos+0.2,pos+i-0.2,color="#80cdc1",lw=2)
		pos+=i
	pos=0
	for label, label_len in zip(labels, label_lens):
		ax.text(-xcax*1.6, pos + label_len / 2, label, color="#80cdc1",ha='right', va='center', rotation=90,
			transform=ax.get_yaxis_transform())
		pos += label_len
	plt.savefig(os.path.join(outdir, "reproducibility.png"), dpi=300, bbox_inches="tight")
	plt.savefig(os.path.join(outdir, "reproducibility.pdf"), dpi=300, bbox_inches="tight")


def binnify(chromsizes_path, binsize, outdir):
	chromsizes = bioframe.read_chromsizes(chromsizes_path, filter_chroms=False)
	bins = bioframe.binnify(chromsizes, binsize)
	bins["mid"] = (bins.start + bins.end) // 2
	outfile = os.path.join(outdir, f"genome_bin{binsize}.bed.gz")
	bins.to_csv(outfile, sep="\t", compression="gzip", index=False, header=False)
	return outfile

def main(myargs):
	#genome_bin = binnify(myargs.chromsize, myargs.resolution, myargs.outdir);exit()

	contact_dict = OrderedDict()
	for line in open(myargs.contact_fofn):
		dataset, contact = line.strip().split()
		contact_dict[dataset] = contact
	
	ds_pairs = list(permutations(contact_dict.keys(), 2))
	#get_reproduce(ds_pairs, contact_dict, genome_bin, myargs.genomedisco, myargs.outdir)
	get_plot(ds_pairs, contact_dict.keys(), myargs.outdir,[i for i in myargs.labels.split(",")],[int(i) for i in myargs.lens.split(",")])
	
def get_myargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-c", "--contact_fofn", required=True, 
			help="contact file of filename, format: <dataset_name>[TAB]<contact_file>. [required]")
	parser.add_argument("-s", "--chromsize", required=True,
			help="chromosome size file, format: <chrom_name>[TAB]<chrom_len>. [required]")
	parser.add_argument("-r", "--resolution", required=True, type=int, 
			help="contact file resolution, used to generate genome-bin file [required]")	
	parser.add_argument("--genomedisco", 
			default="/public/frasergen/PUB/software/python/Python-2.7.9/bin/genomedisco",
			help="genomedisco program path, default: /public/frasergen/PUB/software/python/Python-2.7.9/bin/genomedisco")
	parser.add_argument("--outdir", default="./", 
			help="output dir. [default:./]")
	parser.add_argument("--labels", required=True, 
			help="sample group labels list [required]")
	parser.add_argument("--lens", required=True, 
			help="sample group lengths list [required]")
	return parser.parse_args()

if __name__ == "__main__":
	main(get_myargs())
