import argparse, shutil
import yaml, time
import pandas as pd
from pathlib import Path
from PIL import Image

def resize(infile, outfile, scale=3):
	print(f"[convert] {infile} => {outfile}")
	in_image = Image.open(infile)
	out_width, out_height = in_image.width // scale, in_image.height // scale
	out_image = in_image.resize((out_width, out_height))
	out_image.save(outfile)

def file_operate(inpath, outpath, operate, scale=3):

	if isinstance(inpath, list):
		if not outpath.is_dir(): 
			print(inpath, outpath, operate)
			exit("[ERROR] cannot cp or mv files to one file !")
		for afile in inpath:
			if operate == "cp":
				shutil.copy(afile, outpath)
			elif operate == "resize":
				if afile.suffix.strip(".") not in ["png", "tif", "jpg"]: continue
				afile_base = Path(afile).name
				resize(afile, Path(outpath).joinpath(afile_base), scale)

	elif Path(inpath).is_dir():
		if outpath.exists() and not outpath.is_dir():
			exit("[ERROR] cannot cp or mv dir to one file !")
		if operate == "cp": 
			shutil.copytree(inpath, outpath)

	elif Path(inpath).is_file():
		if Path(outpath).exists() and Path(outpath).is_dir():
			afile_base = Path(inpath).name
			outfile = Path(outpath).joinpath(afile_base)
		else:
			outfile = Path(outpath)
			desdir = Path(outpath).parent
			desdir.mkdir(parents=True, exist_ok=True)

		if operate == "cp":
			shutil.copy(Path(inpath), outfile)
		elif operate == "resize":
			resize(Path(inpath), outfile, scale)

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

	with open(tabdir.joinpath("contract_info.txt"), "w") as outbuff:
		print("合同信息", "合同内容", sep="\t", file=outbuff)
		print("项目编号", conf_dict.get("contract"), sep="\t", file=outbuff)
		print("项目名称", conf_dict.get("title"), sep="\t", file=outbuff)
		print("项目类型", "单样本Hi-C互作标准分析报告", sep="\t", file=outbuff)
		print("报告日期", time.strftime("%Y-%m-%d"), sep="\t", file=outbuff)
	
	samples, datasets = list(), dict()
	for spl_dict in conf_dict.get("samplenames"):
		spl, ds_all = list(spl_dict.items())[0]
		samples.append(spl)
		for ds in ds_all.replace(" ", "").split(","):
			datasets[ds] = spl

	basic_reslu = conf_dict.get("basic_res")
	# step01: QC stat
	prepare_fig01 = figdir.joinpath("01.reads_clean")
	prepare_tab01 = tabdir.joinpath("01.reads_clean")
	result_01 = result_dir.joinpath("01.reads")
	for step1_dir in prepare_fig01, prepare_tab01, result_01:
		step1_dir.mkdir(parents=True, exist_ok=True)
	
	clean_stat = list()
	for ds in datasets:
		raw_clean_dir = indir.joinpath("02.reads", ds, f"{ds}_clean")
		raw_clean_fig = list(raw_clean_dir.glob(f"{ds}.sequencing_quality.p*"))
		file_operate(raw_clean_fig, result_01, "cp")
		file_operate(raw_clean_fig, prepare_fig01, "resize")

		raw_clean_stat = raw_clean_dir.joinpath("quality_stat.txt")
		intable = pd.read_csv(raw_clean_stat, sep="\t")
		intable.rename(columns={"SampleName": "DataSet"}, inplace=True)
		intable.insert(loc=0, column="SampleName", value=datasets[ds])
		clean_stat.append(intable)

	outfile = prepare_tab01.joinpath("sequencing_quality.txt")
	pd.concat(clean_stat).to_csv(outfile, sep="\t", index=False)
	file_operate(outfile, result_01, "cp")

	# step02: alignment
	in_aln_dir = indir.joinpath("03.align")
	prepare_fig02 = figdir.joinpath("02.align")
	prepare_tab02 = tabdir.joinpath("02.align")
	result_02 = result_dir.joinpath("02.align")
	for step2_dir in prepare_fig02, prepare_tab02, result_02:
		step2_dir.mkdir(parents=True, exist_ok=True)

	for ds in datasets:
		ds_stat = in_aln_dir.joinpath(f"{ds}.pair.stat.xls")
		file_operate(ds_stat, result_02, "cp")
		
	aln_stat = pd.DataFrame()
	for spl in samples:
		reslu_assess = list(in_aln_dir.glob(f"{spl}.reslu.p*"))
		file_operate(reslu_assess, result_02, "cp")
		file_operate(reslu_assess, prepare_fig02, "resize")
		
		aln_qc = in_aln_dir.joinpath(f"{spl}_report", "plots")
		proportion = list(aln_qc.glob(f"{spl}.proportion.p*"))
		file_operate(proportion, result_02, "cp")
		file_operate(proportion, prepare_fig02, "resize", 6)
	
		log10prob = list(aln_qc.glob(f"{spl}.log10prob.p*"))
		file_operate(log10prob, result_02, "cp")
		file_operate(log10prob, prepare_fig02, "resize")
		
		raw_aln_stat = in_aln_dir.joinpath(f"{spl}_report", f"{spl}.summary.out")
		intable = pd.read_csv(raw_aln_stat, sep="\t", names=("Category", spl))
		mapped=pd.concat([pd.read_table(in_aln_dir.joinpath(f"{i}.pair.stat"),nrows=19,header=None,index_col=0) for i,j in datasets.items() if j==spl],axis=1)
		map_read=mapped.loc["total"].sum()-mapped.loc["pair_types/NM"].sum()-mapped.loc["pair_types/NN"].sum()
		mapped=pd.DataFrame(["% mapped reads","%.3f" % (map_read/mapped.loc["total"].sum()*100)],index=intable.columns).T
		intable=pd.concat([mapped,intable],axis=0,ignore_index=True)
		if aln_stat.empty: aln_stat = intable
		else: aln_stat = aln_stat.merge(intable, on="Category")

	outfile = prepare_tab02.joinpath("reads_alignment.txt")
	aln_stat = aln_stat[(aln_stat.Category!="convergence") & (aln_stat.Category!="slope")]
	aln_stat.T.to_csv(outfile, sep="\t", header=False)
	file_operate(outfile, result_02, "cp")

	# step03. sample matrix
	in_matrix_dir = indir.joinpath("04.matrix")
	prepare_fig03 = figdir.joinpath("03.heatmap")
	prepare_tab03 = figdir.joinpath("03.heatmap")
	for step3_dir in prepare_fig03, prepare_tab03:
		step3_dir.mkdir(parents=True, exist_ok=True)
	
	genome_res = conf_dict.get("genome_res")
	cpt_res = conf_dict.get("genome_res")
	for spl in samples:
		genome_heatmap = list(in_matrix_dir.glob(f"{spl}_{genome_res}.genome.p*"))
		file_operate(genome_heatmap, result_02, "cp")
		file_operate(genome_heatmap, prepare_fig03, "resize", 5)

		bychr_all_heatmap = in_matrix_dir.joinpath(f"{spl}_{cpt_res}")
		genome_outdir = result_02.joinpath(f"{spl}_{cpt_res}")
		file_operate(bychr_all_heatmap, genome_outdir, "cp")
		if genome_outdir.joinpath("bychr").is_dir():
			shutil.rmtree(genome_outdir.joinpath("bychr"))
	
		chr1_heatmap = in_matrix_dir.joinpath(f"{spl}_{cpt_res}", f"{spl}_{cpt_res}.chr1.png")	
		file_operate(chr1_heatmap, prepare_fig03, "resize", 4)

		spl_cool = in_matrix_dir.joinpath(f"{spl}.basic.cool")
		file_operate(spl_cool, result_02.joinpath(f"{spl}.{basic_reslu}.cool"), "cp")

	reproduce_res = conf_dict.get("reproduce_res")
	reproduce_png = in_matrix_dir.joinpath(f"reproduce_{reproduce_res}", "reproducibility.png")
	file_operate(reproduce_png, result_02, "cp")
	file_operate(reproduce_png, prepare_fig03, "resize")
	reproduce_pdf = in_matrix_dir.joinpath(f"reproduce_{reproduce_res}", "reproducibility.pdf")
	file_operate(reproduce_pdf, result_02, "cp")
	
	# step04. compartment
	cptmt_dir = indir.joinpath("05.compartment")
	cptmt_method = conf_dict.get("compartment").lower()
	cptmt_res = conf_dict.get("compartment_res")
	cptmt_subdir_dict = {"cooltools": "02.PCA", "cscore": "03.Cscore"}
	cptmt_subdir = cptmt_subdir_dict.get(cptmt_method)
	result_04 = result_dir.joinpath("03.compartment")
	prepare_fig04 = figdir.joinpath("04.compartment")
	prepare_tab04 = tabdir.joinpath("04.compartment")
	for step4_dir in result_04, prepare_fig04, prepare_tab04:
		step4_dir.mkdir(parents=True, exist_ok=True)

	gc_list, gene_list, len_list = list(), list(), list()
	for spl in samples:
		bychr_matrix = list(cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
					f"{spl}_{cptmt_res}").glob(f"{spl}_{cptmt_res}.*.p*"))
		file_operate(bychr_matrix, result_04, "cp")
		chr1_matrix = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
					f"{spl}_{cptmt_res}", f"{spl}_{cptmt_res}.chr1.png")
		file_operate(chr1_matrix, prepare_fig04, "resize", 6)

		len_fig = list(cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}").glob( 
					f"{spl}_{cptmt_res}.compartment_length.p*"))
		file_operate(len_fig, result_04, "cp")
		file_operate(len_fig, prepare_fig04, "resize")

		gc_png = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}",
					"Compartment_GC.boxplot.png")
		gc_bname = f"{spl}_{cptmt_res}.compartment.GC.png"
		file_operate(gc_png, result_04.joinpath(gc_bname), "cp")
		file_operate(gc_png, prepare_fig04.joinpath(gc_bname), "resize")
		gc_pdf = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}",
					"Compartment_GC.boxplot.pdf")
		gc_bname = f"{spl}_{cptmt_res}.compartment.GC.pdf"
		file_operate(gc_pdf, result_04.joinpath(gc_bname), "cp")

		cptmt_gene = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
					"Compartment_gene_density.boxplot.png")
		gene_bname = f"{spl}_{cptmt_res}.compartment.gene.png"
		file_operate(cptmt_gene, result_04.joinpath(gene_bname), "cp")
		file_operate(cptmt_gene, prepare_fig04.joinpath(gene_bname), "resize")

		gene_pdf = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
					"Compartment_gene_density.boxplot.pdf")
		gene_bname = f"{spl}_{cptmt_res}.compartment.gene.pdf"
		file_operate(gene_pdf, result_04.joinpath(gene_bname), "cp")
		
		gene_ab = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
								f"{spl}_{cptmt_res}.gc_gene.bed")
		file_operate(gene_ab, result_04, "cp")
		
		cptmt_out = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
								f"{spl}_{cptmt_res}.compartment.xls")
		file_operate(cptmt_out, result_04, "cp")

		gc_cptmt_stat = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
								"Compartment_GC.txt")
		gc_table = pd.read_csv(gc_cptmt_stat, sep="\t")
		gc_table.insert(0, "Sample", spl)
		gc_list.append(gc_table)
		
		gene_cptmt_stat = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
								"Compartment_gene_density.txt")
		gene_table = pd.read_csv(gene_cptmt_stat, sep="\t")
		gene_table.insert(0, "Sample", spl)
		gene_list.append(gene_table)

		cptmt_len_stat = cptmt_dir.joinpath(cptmt_subdir, f"{spl}_{cptmt_res}", 
								f"{spl}_{cptmt_res}.compartment_length.txt")
		len_table = pd.read_csv(cptmt_len_stat, sep="\t")
		len_table.insert(0, "Sample", spl)
		len_list.append(len_table)
	
	cptmt_gc_report = prepare_tab04.joinpath("Compartment_GC.stat.xls")
	pd.concat(gc_list).to_csv(cptmt_gc_report, sep="\t", index=False)
	file_operate(cptmt_gc_report, result_04, "cp")
	cptmt_gene_report = prepare_tab04.joinpath("Compartment_gene_density.stat.xls")
	pd.concat(gene_list).to_csv(cptmt_gene_report, sep="\t", index=False)
	file_operate(cptmt_gene_report, result_04, "cp")
	cptmt_len_report = prepare_tab04.joinpath("Compartment_length.stat.xls")
	pd.concat(len_list).to_csv(cptmt_len_report, sep="\t", index=False)
	file_operate(cptmt_len_report, result_04, "cp")

	# step05. TAD
	tad_dir = indir.joinpath("06.tad")
	tad_method = conf_dict.get("tad").lower()
	tad_res = conf_dict.get("tad_res")
	tad_subdir_dict = {"cooltools": "02.cooltools", "hicfindtads": "03.hicFindTAD", 
				"cworld": "04.cworld"}
	tad_subdir = tad_subdir_dict.get(tad_method)
	result_05 = result_dir.joinpath("04.TAD")
	prepare_fig05 = figdir.joinpath("05.TAD")
	prepare_tab05 = tabdir.joinpath("05.TAD")
	for step5_dir in result_05, prepare_fig05, prepare_tab05:
		step5_dir.mkdir(parents=True, exist_ok=True)
	tad_gc_list, tad_gene_list, tad_len_list = list(), list(), list()
	for spl in samples:
		tad_bed = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", "TAD_genome.bed")
		file_operate(tad_bed, result_05.joinpath(f"{spl}_{tad_res}.TAD.bed"), "cp")
		
		tad_gc_png = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", 
									"TAD_border_and_inter_GC.boxplot.png")
		file_operate(tad_gc_png, result_05.joinpath(f"{spl}_{tad_res}.TAD_GC.png"), "cp")
		file_operate(tad_gc_png, prepare_fig05.joinpath(f"{spl}_{tad_res}.TAD_GC.png"), "resize")
		tad_gc_pdf = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", 
									"TAD_border_and_inter_GC.boxplot.pdf")
		file_operate(tad_gc_pdf, result_05.joinpath(f"{spl}_{tad_res}.TAD_GC.pdf"), "cp")
		
		tad_gene_png = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", 
									"TAD_border_and_inter_gene_density.boxplot.png")
		file_operate(tad_gene_png, result_05.joinpath(f"{spl}_{tad_res}.TAD_Gene.png"), "cp")
		file_operate(tad_gene_png, prepare_fig05.joinpath(f"{spl}_{tad_res}.TAD_Gene.png"), "resize")
		tad_gene_pdf = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", 
									"TAD_border_and_inter_gene_density.boxplot.pdf")
		file_operate(tad_gene_pdf, prepare_fig05.joinpath(f"{spl}_{tad_res}.TAD_Gene.pdf"), "cp")
		
		out_display = prepare_fig05.joinpath("Display")
		out_display.mkdir(parents=True, exist_ok=True)
		display_png = list(tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", 
									"TAD_Display").glob("*.png"))
		file_operate(display_png, out_display, "resize", 8)

		boundary_motif = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", "TAD_boundary_motif.txt")
		bmotif_dist = result_05.joinpath(f"{spl}_{tad_res}.boundary_motif.txt")
		file_operate(boundary_motif, bmotif_dist, "cp")
		boundary_motif_png = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", f"{spl}.{tad_res}_boundary_motif5.png")
		file_operate(boundary_motif_png, prepare_fig05, "cp")
		file_operate(boundary_motif_png, result_05, "cp")
		boundary_motif_pdf = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", f"{spl}.{tad_res}_boundary_motif5.pdf")
		file_operate(boundary_motif_pdf, result_05, "cp")

		gc_tad_stat = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", "TAD_border_and_inter_GC.txt")
		gc_tad_table = pd.read_csv(gc_tad_stat, sep="\t")
		gc_tad_table.insert(0, "Sample", spl)
		tad_gc_list.append(gc_tad_table)

		gene_tad_stat = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", "TAD_border_and_inter_gene_density.txt")
		gene_tad_table = pd.read_csv(gene_tad_stat, sep="\t")
		gene_tad_table.insert(0, "Sample", spl)
		tad_gene_list.append(gene_tad_table)

		len_tad_stat = tad_dir.joinpath(tad_subdir, f"{spl}.{tad_res}", "Domain_length.stat.txt")
		len_tad_table = pd.read_csv(len_tad_stat, sep="\t")
		len_tad_table.insert(0, "Sample", spl)
		tad_len_list.append(len_tad_table)

	tad_gc_report = prepare_tab05.joinpath("TAD_GC.stat.xls")
	pd.concat(tad_gc_list).to_csv(tad_gc_report, sep="\t", index=False)
	file_operate(tad_gc_report, result_05, "cp")

	tad_gene_report = prepare_tab05.joinpath("TAD_gene_density.stat.xls")
	pd.concat(tad_gene_list).to_csv(tad_gene_report, sep="\t", index=False)
	file_operate(tad_gene_report, result_05, "cp")

	tad_len_report = prepare_tab05.joinpath("TAD_length.stat.xls")
	pd.concat(tad_len_list).to_csv(tad_len_report, sep="\t", index=False)
	file_operate(tad_len_report, result_05, "cp")
	
	# step06. loop
	result_06 = result_dir.joinpath("05.loop")
	prepare_fig06 = figdir.joinpath("06.loop")
	prepare_tab06 = tabdir.joinpath("06.loop")
	loop_stat = prepare_tab06.joinpath("loop_stat.txt")
	for step6_dir in result_06, prepare_fig06, prepare_tab06:
		step6_dir.mkdir(parents=True, exist_ok=True)

	loop_method = conf_dict.get("loop").lower()
	loop_res = str(conf_dict.get("loop_res")).replace(" ", "").replace(",", "_")
	
	with open(loop_stat, "w") as outbuff:
		if loop_method == "fithic":
			if "_" in loop_res: exit("\033[91m\n[ERROR] fithic-based report support only one resolution result \033[0m")
			print("Sample", "Cis_loops", "Trans_loops", sep="\t", file=outbuff)
			loop_subdir = indir.joinpath("07.loops", "03.fithic")
			for spl in samples:
				cis_png = loop_subdir.joinpath(f"{spl}.{loop_res}", "04.circos", "cis_interaction.png")
				file_operate(cis_png, result_06.joinpath(f"{spl}_{loop_res}.cis_interaction.png"), "cp")
				file_operate(cis_png, prepare_fig06.joinpath(f"{spl}_{loop_res}.cis_interaction.png"), "resize", 6)
				cis_svg = loop_subdir.joinpath(f"{spl}.{loop_res}", "04.circos", "cis_interaction.svg")
				file_operate(cis_svg, result_06.joinpath(f"{spl}_{loop_res}.cis_interaction.svg"), "cp")
			
				trans_png = loop_subdir.joinpath(f"{spl}.{loop_res}", "04.circos", "trans_interaction.png")
				file_operate(trans_png, result_06.joinpath(f"{spl}_{loop_res}.trans_interaction.png"), "cp")
				file_operate(trans_png, prepare_fig06.joinpath(f"{spl}_{loop_res}.trans_interaction.png"), "resize", 6)
				trans_svg = loop_subdir.joinpath(f"{spl}.{loop_res}", "04.circos", "trans_interaction.svg")
				file_operate(trans_svg, result_06.joinpath(f"{spl}_{loop_res}.trans_interaction.svg"), "cp")

				cis_bedpe = loop_subdir.joinpath(f"{spl}.{loop_res}", "04.circos", f"{spl}_{loop_res}.cis_interact.bedpe")
				file_operate(cis_bedpe, result_06, "cp")
				cis_annot = loop_subdir.joinpath(f"{spl}.{loop_res}", "04.circos", f"{spl}_2000.cis_annot.xls")
				file_operate(cis_annot, result_06, "cp")
				trans_bedpe = loop_subdir.joinpath(f"{spl}.{loop_res}", "04.circos", f"{spl}_{loop_res}.trans_interact.bedpe")
				file_operate(trans_bedpe, result_06, "cp")
				trans_annot = loop_subdir.joinpath(f"{spl}.{loop_res}", "04.circos", f"{spl}_2000.trans_annot.xls")
				file_operate(trans_annot, result_06, "cp")
				cis_loopnum = sum([1 for i in open(cis_bedpe) if not i.startswith("#")])
				trans_loopnum = sum([1 for i in open(trans_bedpe) if not i.startswith("#")])
				print(spl, cis_loopnum, trans_loopnum, sep="\t", file=outbuff)

				

		if loop_method == "mustache":
			print("Sample", "Cis_loops", sep="\t", file=outbuff)
			loop_subdir = indir.joinpath("07.loops", "02.mustache")
			for spl in samples:
				cis_png = loop_subdir.joinpath(spl, f"{spl}_circos", "cis_interaction.png")
				file_operate(cis_png, result_06.joinpath(f"{spl}.cis_interaction.png"), "cp")
				file_operate(cis_png, prepare_fig06.joinpath(f"{spl}.cis_interaction.png"), "resize", 6)
				cis_svg = loop_subdir.joinpath(spl, f"{spl}_circos", "cis_interaction.svg")
				file_operate(cis_svg, result_06.joinpath(f"{spl}.cis_interaction.svg"), "cp")
				
				cis_loops = loop_subdir.joinpath(spl, f"{spl}.combine.bedpe")
				file_operate(cis_loops, result_06, "cp")
				loop_num = sum([1 for i in open(cis_loops) if not i.startswith("#")])
				print(spl, loop_num, sep="\t",file=outbuff)

				loop_annot = loop_subdir.joinpath(spl, f"{spl}_annot.xls")
				file_operate(loop_annot, result_06, "cp")
	file_operate(loop_stat, result_06, "cp")
		
def get_myargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--indir", required=True, 
			help="input single sample workdir [required]")
	parser.add_argument("-y", "--yaml_file", required=True, 
			help="single sample pipeline yaml config file [required]")
	parser.add_argument("-o", "--outdir", default="./", 
			help="output dir [default: ./")
	return parser.parse_args()


if __name__ == "__main__":
	main(get_myargs())	

