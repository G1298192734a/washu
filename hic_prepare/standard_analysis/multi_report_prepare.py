import argparse, shutil
import yaml, time, os
import pandas as pd
from pathlib import Path
from single_report_prepare import file_operate


def main(myargs):
	indir = Path(myargs.indir)
	outdir = Path(myargs.outdir)
	prepare_dir = outdir.joinpath("prepare_report")
	if prepare_dir.exists(): shutil.rmtree(prepare_dir)
	prepare_dir.mkdir(parents=True, exist_ok=True)
	
	result_dir = outdir.joinpath("result")
	if result_dir.exists(): shutil.rmtree(result_dir)
	result_dir.mkdir(parents=True, exist_ok=True)
	
	figdir = prepare_dir.joinpath("figure")
	tabdir = prepare_dir.joinpath("table")
	tabdir.mkdir(parents=True, exist_ok=True)

	conf_dict = yaml.safe_load(open(myargs.yaml_file))
	cptmt_res = conf_dict.get("compartment_res")
	with open(tabdir.joinpath("contract_info.txt"), "w") as outbuff:
		print("合同信息", "合同内容", sep="\t", file=outbuff)
		print("项目编号", conf_dict.get("contract"), sep="\t", file=outbuff)
		print("项目名称", conf_dict.get("title"), sep="\t", file=outbuff)
		print("项目类型", "多样本Hi-C互作标准分析报告", sep="\t", file=outbuff)
		print("报告日期", time.strftime("%Y-%m-%d"), sep="\t", file=outbuff)

	spl_compare = conf_dict.get("sample_minus").strip().split(",")
	with open(tabdir.joinpath("sample_info.txt"), "w") as outbuff:
		print("sample_pairs", "control", "treat", sep="\t", file=outbuff)
		for spl_pair in spl_compare:
			treat, control = spl_pair.split(".minus.")
			print(spl_pair, control, treat, sep="\t", file=outbuff)

	
	# step01. differential heatmap
	prepare_fig01 = figdir.joinpath("01.diff_heatmap")
	result_01 = result_dir.joinpath("01.diff_heatmap")
	for step1_dir in prepare_fig01, result_01:
		step1_dir.mkdir(parents=True, exist_ok=True)

	for spl_pair in spl_compare:
		diff_genome = indir.joinpath("01.diff_heatmap", "01.genome_diffhm", spl_pair, f"{spl_pair}.genome.png")
		file_operate(diff_genome, prepare_fig01, "resize")
		diff_chr1 = indir.joinpath("01.diff_heatmap", "02.bychr_diffhm", spl_pair, f"{spl_pair}.chr1.png")
		file_operate(diff_chr1, prepare_fig01, "resize")
		decay_curve = indir.joinpath("01.diff_heatmap", "03.diff_pscurve", "decay_curve.png")
		file_operate(decay_curve, prepare_fig01, "resize")

		result01_subdir = result_01.joinpath(spl_pair)
		result01_subdir.mkdir(parents=True, exist_ok=True)
		gdhms = list(indir.joinpath("01.diff_heatmap", "01.genome_diffhm", spl_pair).glob(f"{spl_pair}.*"))
		file_operate(gdhms, result01_subdir, "cp")
		
		bychr_fs = list(indir.joinpath("01.diff_heatmap", "02.bychr_diffhm", spl_pair).glob("*"))
		file_operate(bychr_fs, result01_subdir, "cp")
		pscurve = list(indir.joinpath("01.diff_heatmap", "03.diff_pscurve").glob("*"))
		file_operate(pscurve, result01_subdir, "cp")

	# step02. differential compartment
	prepare_fig02 = figdir.joinpath("02.diff_compartment")
	prepare_tab02 = tabdir.joinpath("02.diff_compartment")
	result_02 = result_dir.joinpath("02.diff_compartment")
	for step2_dir in prepare_fig02, prepare_tab02, result_02:
		step2_dir.mkdir(parents=True, exist_ok=True)
	
	switch_table = list()
	for spl_pair in spl_compare:
		cptmt_figs = list(indir.joinpath("02.diff_compartment", spl_pair).glob(f"*.png"))
		file_operate(cptmt_figs, prepare_fig02, "resize")
		
		switch_plot = indir.joinpath("02.diff_compartment", spl_pair, f"switch_plot.{cptmt_res}", f"{spl_pair}.chr1.png")
		file_operate(switch_plot, prepare_fig02, "resize")
		A2B_enrich_figs = list(indir.joinpath("02.diff_compartment", spl_pair, f"{spl_pair}.A2B").glob("*.png"))
		file_operate(A2B_enrich_figs, prepare_fig02, "cp")
		#B2A_enrich_figs = list(indir.joinpath("02.diff_compartment", spl_pair, f"{spl_pair}.B2A").glob("*.png"))
		#file_operate(B2A_enrich_figs, prepare_fig02, "resize")
		cptmt_all = list(indir.joinpath("02.diff_compartment", "All").glob("*.png"))
		file_operate(cptmt_all, prepare_fig02, "resize")
		
		switch_file = indir.joinpath("02.diff_compartment", spl_pair, f"{spl_pair}.switch_stat.xls")
		intable = pd.read_csv(switch_file, sep="\t")
		switch_table.append(intable)
		
		splpair_dir = indir.joinpath("02.diff_compartment", spl_pair)
		file_operate(splpair_dir, result_02.joinpath(spl_pair), "cp")

		splall = list(indir.joinpath("02.diff_compartment", "All").glob("*"))
		file_operate(splall, result_02, "cp")

	switch_alltab = prepare_tab02.joinpath("switch_stat.xls")
	pd.concat(switch_table).to_csv(switch_alltab, sep="\t", index=False)
	
	# step03. differential TAD
	prepare_fig03 = figdir.joinpath("03.diff_TAD")
	prepare_tab03 = tabdir.joinpath("03.diff_TAD")
	result_03 = result_dir.joinpath("03.diff_TAD")
	for step3_dir in prepare_fig03, prepare_tab03, result_03:
		step3_dir.mkdir(parents=True, exist_ok=True)
	
	classify_table = list()
	for spl_pair in spl_compare:
		tad_figs = list(indir.joinpath("03.diff_TAD", "01.cooltools", spl_pair).glob("*.png"))
		file_operate(tad_figs, prepare_fig03, "resize")
		specific_enrich = list(indir.joinpath("03.diff_TAD", "01.cooltools", spl_pair).glob("*/*.png"))
		file_operate(specific_enrich, prepare_fig03, "cp")
		
	
		classify = indir.joinpath("03.diff_TAD", "01.cooltools", spl_pair, f"{spl_pair}.classify.stat.xls")
		intable = pd.read_csv(classify, sep="\t", index_col=0).T
		intable = intable.reset_index().rename(columns={"index": "Sample"})
		intable.insert(0, "Sample_compare", spl_pair)
		classify_table.append(intable)

		splpair_dir = indir.joinpath("03.diff_TAD", "01.cooltools", spl_pair)
		file_operate(splpair_dir, result_03.joinpath(spl_pair), "cp")

	classify_alltab = prepare_tab03.joinpath("classify_stat.xls")
	pd.concat(classify_table).to_csv(classify_alltab, sep="\t", index=False)
		
	# step04. differential loops
	prepare_fig04 = figdir.joinpath("04.diff_loop")
	prepare_tab04 = tabdir.joinpath("04.diff_loop")
	result_04 = result_dir.joinpath("04.diff_loop")
	for step4_dir in prepare_fig04, prepare_tab04, result_04:
		step4_dir.mkdir(parents=True, exist_ok=True)

	loop_caller = conf_dict.get("loop")
	if loop_caller == "mustache":
		looplen_allpng = indir.joinpath("04.diff_loop", "01.mustache", "All", "loop_length_all.png")
		file_operate(looplen_allpng, prepare_fig04, "resize")
		file_operate(looplen_allpng, result_04, "cp")
		looplen_allpdf = indir.joinpath("04.diff_loop", "01.mustache", "All", "loop_length_all.pdf")
		file_operate(looplen_allpdf, result_04, "cp")

		diffloop_table = list()
		for spl_pair in spl_compare:
			apa_fig = indir.joinpath("04.diff_loop", "01.mustache", spl_pair, f"{spl_pair}.APA.png")
			file_operate(apa_fig, prepare_fig04, "resize")
			loopstat = indir.joinpath("04.diff_loop", "01.mustache", spl_pair, f"{spl_pair}_loop_stat.xls")
			intable = pd.read_csv(loopstat, sep="\t")
			intable.insert(0, "Sample_compare", spl_pair)
			diffloop_table.append(intable)
	
			treat_specific = spl_pair + "-" + spl_pair.split(".minus.")[0] + "_specific"
			treat_enrich = list(indir.joinpath("04.diff_loop", "01.mustache", 
								spl_pair, treat_specific).glob("*.png"))
			file_operate(treat_enrich, prepare_fig04, "cp")
	
			control_specific = spl_pair + "-" + spl_pair.split(".minus.")[1] + "_specific"
			control_enrich = list(indir.joinpath("04.diff_loop", "01.mustache", 
								spl_pair, control_specific).glob("*.png"))
			file_operate(control_enrich, prepare_fig04, "cp")
			
			splpair_dir = indir.joinpath("04.diff_loop", "01.mustache", spl_pair)
			file_operate(splpair_dir, result_04.joinpath(spl_pair), "cp")
		diffloop_alltab = prepare_tab04.joinpath("diff_loop_stat.xls")
		pd.concat(diffloop_table).to_csv(diffloop_alltab, sep="\t", index=False)

	elif loop_caller == "fithic":
		looplen_allpng = indir.joinpath("04.diff_loop", "02.fithic", "All", "loop_length_all.png")
		file_operate(looplen_allpng, prepare_fig04, "resize")
		file_operate(looplen_allpng, result_04, "cp")
		looplen_allpdf = indir.joinpath("04.diff_loop", "02.fithic", "All", "loop_length_all.pdf")
		file_operate(looplen_allpdf, result_04, "cp")
		
		diffloop_alltab = prepare_tab04.joinpath("diff_loop_stat.xls")
		outbuff = open(diffloop_alltab, "w")
		print("sample_compare", "loops", "cis_or_trans", "loop_number", sep="\t", file=outbuff)
		for spl_pair in spl_compare:
			for cistrans in "cis", "trans":
				for splid in spl_pair.split(".minus."):
					fithic_bedpe = indir.joinpath("04.diff_loop", "02.fithic", spl_pair, 
							f"{spl_pair}-{splid}_specific_{cistrans}.bedpe")
					loop_num = sum(1 for line in open(fithic_bedpe))
					print(spl_pair, f"{splid}_specific", cistrans, f"{loop_num:,}", sep="\t", file=outbuff)
					file_operate(fithic_bedpe, result_04, "cp")	
					
					genelist = indir.joinpath("04.diff_loop", "02.fithic", spl_pair,
								f"{spl_pair}-{splid}_specific_{cistrans}.genelist")
					file_operate(genelist, result_04, "cp")

					go_kegg = indir.joinpath("04.diff_loop", "02.fithic", spl_pair, 
								f"{spl_pair}-{splid}_specific_{cistrans}")
					file_operate(go_kegg, result_04.joinpath(f"{spl_pair}-{splid}_specific_{cistrans}"), "cp")
				
				share_bedpe = indir.joinpath("04.diff_loop", "02.fithic", spl_pair,
								f"{spl_pair}-share_{cistrans}.bedpe")
				loop_num = sum(1 for line in open(share_bedpe))
				print(spl_pair, "share", cistrans, f"{loop_num:,}", sep="\t", file=outbuff)
				file_operate(share_bedpe, result_04, "cp")
			
				genelist = indir.joinpath("04.diff_loop", "02.fithic", spl_pair,
								f"{spl_pair}-share_{cistrans}.genelist")
				file_operate(genelist, result_04, "cp")

				go_kegg = indir.joinpath("04.diff_loop", "02.fithic", spl_pair, 
								f"{spl_pair}-share_{cistrans}")
				file_operate(go_kegg, result_04.joinpath(f"{spl_pair}-share_{cistrans}"), "cp")
				
		outbuff.close()		
		cis_share = list(indir.joinpath("04.diff_loop", "02.fithic", spl_pair,
								f"{spl_pair}-share_cis").glob("*.png"))
		file_operate(cis_share, prepare_fig04, "cp")

		for cis_trans in "cis", "trans":
			venn_png = indir.joinpath("04.diff_loop", "02.fithic", spl_pair, f"{spl_pair}-{cis_trans}.venn.png")
			file_operate(venn_png, prepare_fig04, "resize")
			file_operate(venn_png, result_04, "cp")	
			venn_pdf = indir.joinpath("04.diff_loop", "02.fithic", spl_pair, f"{spl_pair}-{cis_trans}.venn.pdf")
			file_operate(venn_pdf, result_04, "cp")

	for od in prepare_dir, result_dir:
		for root, dirs, files in os.walk(od):
			for f in files:
				if ".HeatWord." in f:
					os.remove(os.path.join(root, f))

def get_myargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--indir", required=True, 
			help="input multi sample workdir [required]")
	parser.add_argument("-y", "--yaml_file", required=True, 
			help="hic3_biorep pipeline yaml config file [required]")
	parser.add_argument("-o", "--outdir", default="./", 
			help="output dir [default: ./]")
	return parser.parse_args()
	
if __name__ == "__main__":
	main(get_myargs())
