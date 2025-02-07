import sys, snakemake, yaml
import os, time, pysam, math
from os.path import dirname
from PIL import Image
from diff_matrix import parse_genome_cool,parse_chrom_cool
absjoin = lambda *intuple: os.path.abspath(os.path.join(*intuple))

snake_dir = "/public/frasergen/3D/pipeline/Interactome/hic3_workflow"
inyaml = absjoin(snake_dir, "software.yaml")
dbyaml = absjoin(snake_dir, "database.yaml")
scripts_dir = absjoin(snake_dir, "scripts")
src_dir = absjoin(snake_dir, "src")
software = yaml.safe_load(open(inyaml))
database = yaml.safe_load(open(dbyaml))
species_dict = {"a": ("animal_kegg", "animal_nr"), 
				"p": ("plant_kegg", "plant_nr"), 
				"f": ("fungi_kegg", "fungi_nr"), 
				"b": ("bacteria_kegg", "bacteria_nr")}
kegg, nr = species_dict[config.get("species")[0].lower()]
kegg_db, nr_db = database[kegg], database[nr]

samples = config.get("samplenames").replace(" ", "").split(",")
datasets, ds_dict, reads_dict,lens = list(), dict(), dict(),list()
for spl in samples:
	lens.append(len(config[spl]))
	for dataset, pereads in config[spl].items():
		datasets.append(dataset)
		ds_dict[dataset] = spl
		reads_dict[dataset] = pereads

get_reslu = lambda x: str(config.get("resolution").get(x)).replace(" ", "").split(",")
lowest_tad_reslu = max(map(int, get_reslu("tad")))

# make dirs
snakemake.utils.makedirs(["01.ref", "02.reads", "03.align", "04.matrix", "05.compartment", 
				"06.tad", "logdir", absjoin("01.ref", "KEGG"), absjoin("01.ref", "NR"), 
				absjoin("01.ref", "motif_scan")])
sub_dirs = list()
for ds in datasets:
	sub_dir1 = absjoin("02.reads", ds)
	sub_dir2 = absjoin("02.reads", ds, f"{ds}_clean")
	sub_dir3 = absjoin("02.reads", ds, f"{ds}_clean_fastqc")
	sub_dir4 = absjoin("05.compartment", "01.prepare")
	sub_dirs.extend([sub_dir1, sub_dir2, sub_dir3, sub_dir4])

for spl in samples:
	for cptmt_reslu in get_reslu("compartment"):
		pca_dir = absjoin("05.compartment", "02.PCA", f"{spl}_{cptmt_reslu}")
		sub_dirs.append(pca_dir)
		cscore_dir = absjoin("05.compartment", "03.Cscore", f"{spl}_{cptmt_reslu}")
		sub_dirs.append(cscore_dir)
	
	for tad_reslu in get_reslu("tad"):
		findtad = absjoin("06.tad", "02.cooltools", f"{spl}.{tad_reslu}")
		sub_dirs.append(findtad)
		is_dir = absjoin("06.tad", "03.hicFindTAD", f"{spl}.{tad_reslu}")
		sub_dirs.append(is_dir)
		diamond_dir = absjoin("06.tad", "04.cworld", f"{spl}.{tad_reslu}")
		sub_dirs.append(diamond_dir)

snakemake.utils.makedirs(sub_dirs)

# pairqc max logdist
def get_maxlog_dist(inref):
	max_len = max(pysam.FastaFile(inref).lengths)
	return round(math.log10(max_len), 2)

def get_adapter(platform):
	if platform.lower() == "bgi":
		return absjoin(src_dir, "BGI-SEQ-PE.fa")
	elif platform.lower() == "illumina":
		return absjoin(src_dir, "TruSeq3-PE-2.fa")
	else:
		sys.exit("[ERROR] ADAPTER must be either illumina or BGI !")

def reslu_contact(reslu):
	contacts = list()
	for dataset in datasets:
		contact_path = absjoin("04.matrix", f"reproduce_{reslu}", dataset, 
				f"refined.frag.{reslu}.contact.gz")
		contacts.append(contact_path)
	return contacts

rule all:
	input:
		absjoin("01.ref", "KEGG", "genome.ACC2KEGG"),
		absjoin("01.ref", "NR", "genome.ACC2GO"),
		expand(absjoin("02.reads", "{dataset}", "{dataset}_clean", 
				"{dataset}.sequencing_quality.png"), dataset=datasets),
		expand(absjoin("03.align", "{sample}_report", "plots", "{sample}.proportion.png"), sample=samples),
		expand(absjoin("03.align", "{dataset}.nodup.pair.gz"), dataset=datasets),
		expand(absjoin("04.matrix", "reproduce_{reslu}", "work.sh"), reslu=get_reslu("reproduce")),
		#### for sample pair & cool matrix
		expand(absjoin("03.align", "{sample}.merged.pair.gz"), sample=samples),
		expand(absjoin("03.align", "{sample}.reslu.png"), sample=samples),
		expand(absjoin("03.align", "{sample}_report", "{sample}.cis_to_trans.out"), sample=samples),
		expand(absjoin("04.matrix", "{sample}_{genome_reslu}.genome.png"),
						sample=samples, genome_reslu=get_reslu("genome")),
		expand(absjoin("04.matrix", "{sample}_{reslu}", "{sample}.zScore.matrix.gz"), reslu=get_reslu("reproduce"),
						sample=samples), #genome_reslu=get_reslu("genome")),
		#### for sample compartment
		expand(absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", 
						"{sample}_{cptmt_reslu}.bw"), sample=samples, cptmt_reslu=get_reslu("compartment")),
		expand(absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", 
						"{sample}_{cptmt_reslu}.compartment.xls"), sample=samples, 
						cptmt_reslu=get_reslu("compartment")),
		expand(absjoin("05.compartment", "03.Cscore", "{sample}_{cptmt_reslu}", 
							"{sample}_{cptmt_reslu}.bw"), sample=samples, cptmt_reslu=get_reslu("compartment")),
		expand(absjoin("04.matrix","{sample}_{cptmt_reslu}","bychr","{sample}.chr1.zScore.matrix.gz"),sample=samples, cptmt_reslu=get_reslu("compartment")),
		#### for sample TAD
		expand(absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}", "TAD_genome.bed"),
						sample=samples, tad_reslu=get_reslu("tad")),
		expand(absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}", "TAD_genome.bed"),
						sample=samples, tad_reslu=get_reslu("tad")),
		expand(absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}", "TAD_genome.bed"),
						sample=samples, tad_reslu=get_reslu("tad")),
		#### for sample loops
		expand(absjoin("07.loops", "02.mustache", "{sample}", "{sample}.combine.mustache"), 
						sample=samples),
		expand(absjoin("07.loops", "02.mustache", "{sample}", "{sample}_circos", 
						"cis_interaction.png"), sample=samples),
		expand(absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
						"cis_interaction.png"), sample=samples, loop_reslu=get_reslu("loop")),

rule fithic_circos:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		faidx = absjoin("01.ref", "genome.fa.fai"),
		genes = absjoin("01.ref", "genome_gene.bed"),
		cis_fithic = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "02.fithic_cis",
					"{sample}.{loop_reslu}.cis.spline_pass2.res{loop_reslu}.significances.txt.gz"),
		trans_fithic = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "03.fithic_trans",
					"{sample}.{loop_reslu}.trans.spline_pass1.res{loop_reslu}.significances.txt.gz"),
	output:
		cis_interact = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
					"{sample}_{loop_reslu}.cis.interact.gz"),
		cis_bedpe = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos", 
					"{sample}_{loop_reslu}.cis_interact.bedpe"),
		trans_interact = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
					"{sample}_{loop_reslu}.trans.interact.gz"),
		trans_bedpe = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos", 
					"{sample}_{loop_reslu}.trans_interact.bedpe"),
		karyotype = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
						"{sample}_{loop_reslu}.karyotype.txt"),
		gene_density = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
						"{sample}_{loop_reslu}.gene_density.bedgraph"),
		cis_density = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
						"{sample}_{loop_reslu}.cis_density.bedgraph"),
		trans_density = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
						"{sample}_{loop_reslu}.trans_density.bedgraph"),
		cis_links = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
						"{sample}_{loop_reslu}.cis.links"),
		trans_links = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
						"{sample}_{loop_reslu}.trans.links"),
		cis_figure = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos",
						"cis_interaction.png"),
	params:
		pigz = "/public/frasergen/PUB/software/Anaconda/anaconda3-3d/bin/pigz",
		python = software.get("python3"),
		bedtools = software.get("bedtools"),
		#loop_bridge = absjoin(scripts_dir, "loop_bridge.py"),
		loop_bridge="/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/loop_bridge.py",
		fithic_links = absjoin(scripts_dir, "fithic_links.py"),
		run_circos = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/run_circos.py",
		outdir = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos"),
		outpfix = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "04.circos", "{sample}_{loop_reslu}"),
	shell:
		"{params.pigz} -p 30 -cd {input.cis_fithic} |"
		"    awk '$1 == $3 && $5 > 2 && $6 < 0.01 && $7 < 0.01' |"
		"    sort -gk6 | {params.pigz} -p 30 > {output.cis_interact} ;\n"

		"{params.pigz} -p 30 -cd {output.cis_interact} |"
		"""    awk -v OFS="\\t" '{{print "chr"$1,$2-half_res,$2+half_res,"chr"$3,$4-half_res,$4+half_res}}' """
		"""    half_res=$[{wildcards.loop_reslu}/2] |"""
		"      {params.bedtools} sort -faidx {input.faidx} -i stdin > {output.cis_bedpe} ;\n"""

		"{params.python} {params.loop_bridge} --loops_bedpe {output.cis_bedpe}"
		"    --gene_annot {input.genes} --annot_type bed --outpfix {params.outpfix}.cis_annot ;\n"

		"{params.pigz} -p 30 -cd {input.trans_fithic} |"
		"    awk '$1 != $3 && $5 > 2 && $6 < 0.01 && $7 < 0.01' |"
		"    sort -nrk5 | {params.pigz} -p 30 > {output.trans_interact} ;\n"
		
		"{params.pigz} -p 30 -cd {output.trans_interact} |"
		"""    awk -v OFS="\\t" '{{print "chr"$1,$2-half_res,$2+half_res,"chr"$3,$4-half_res,$4+half_res}}' """
		"""    half_res=$[{wildcards.loop_reslu}/2] |"""
		"      {params.bedtools} sort -faidx {input.faidx} -i stdin > {output.trans_bedpe} ;\n"""

		"{params.python} {params.loop_bridge} --loops_bedpe {output.trans_bedpe}"
		"     --gene_annot {input.genes} --annot_type bed --outpfix {params.outpfix}.trans_annot ;\n"

		"{params.python} {params.fithic_links} --chromsize {input.reflen}"
		"    --fithic_cis {output.cis_interact} --fithic_trans {output.trans_interact}"
		"    --resolution {wildcards.loop_reslu} --gene_annot {input.genes}"
		"    --outpfix {params.outpfix} ;\n"
		"export PATH=/work/frasergen/PUB/software/perl/perl-5.34.0/bin:$PATH ;\n"
		"{params.python} {params.run_circos} --method cis --karyotype {output.karyotype}" 
		"    --gene_density {output.gene_density} --link_density {output.cis_density}"
		"    --color_links {output.cis_links} --outdir {params.outdir} ;\n"

		"{params.python} {params.run_circos} --method trans --karyotype {output.karyotype}"
		"    --gene_density {output.gene_density} --link_density {output.trans_density}"
		"    --color_links {output.trans_links} --outdir {params.outdir} ;\n"

rule fithic_loops:
	input:
		loop_cool = absjoin("07.loops", "01.prepare", "{sample}_{loop_reslu}.cool"),
	output:
		contact = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "01.fithic_input",
					"refined.frag.{loop_reslu}.contact.gz"),
		bias = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "01.fithic_input",
					"iced.{loop_reslu}.bias.gz"),
		fragment = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "01.fithic_input",
					"refined.frag.{loop_reslu}.fragment.gz"),
		cis_fithic = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "02.fithic_cis", 
					"{sample}.{loop_reslu}.cis.spline_pass2.res{loop_reslu}.significances.txt.gz"),
		trans_fithic = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "03.fithic_trans", 
					"{sample}.{loop_reslu}.trans.spline_pass1.res{loop_reslu}.significances.txt.gz"),
	params:
		python = software.get("python3"),
		fithic = software.get("fithic"),
		lowerbound = lambda wildcards: 2 * int(wildcards.loop_reslu),
		cooler_fithic = absjoin(scripts_dir, "cooler_fithic2.py"),
		outdir = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "01.fithic_input"),
		cisdir = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "02.fithic_cis"),
		transdir = absjoin("07.loops", "03.fithic", "{sample}.{loop_reslu}", "03.fithic_trans"),
	shell:
		"{params.python} {params.cooler_fithic} --incool {input.loop_cool}"
		"    --outdir {params.outdir} ;\n"

		"{params.python} {params.fithic} --interactions {output.contact}"
		"    --fragments {output.fragment} --biases {output.bias}"
		"    --lowerbound {params.lowerbound} --outdir {params.cisdir}"
		"    --resolution {wildcards.loop_reslu} --visual --passes 2"
		"    --contactType intraOnly --lib {wildcards.sample}.{wildcards.loop_reslu}.cis ;\n"

		"{params.python} {params.fithic} --interactions {output.contact}"
		"    --fragments {output.fragment} --outdir {params.transdir}"
		"    --resolution {wildcards.loop_reslu} --visual --contactType interOnly"
		"    --lib {wildcards.sample}.{wildcards.loop_reslu}.trans ;\n"		

rule mustache_circos:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		genes = absjoin("01.ref", "genome_gene.bed"),
		loop_combine = absjoin("07.loops", "02.mustache", "{sample}", 
							"{sample}.combine.mustache"),
	output: 
		karyotype = absjoin("07.loops", "02.mustache", "{sample}", 
							"{sample}_circos", "{sample}.karyotype.txt"),
		gene_density = absjoin("07.loops", "02.mustache", "{sample}", 
							"{sample}_circos", "{sample}.gene_density.bedgraph"),
		link_density = absjoin("07.loops", "02.mustache", "{sample}", 
							"{sample}_circos", "{sample}.cis_density.bedgraph"),
		color_links = absjoin("07.loops", "02.mustache", "{sample}", 
							"{sample}_circos", "{sample}.cis.links"),
		circos_figure = absjoin("07.loops", "02.mustache", "{sample}", 
							"{sample}_circos", "cis_interaction.png"),
	params:
		python3 = software.get("python3"),
		mustache_circos = absjoin(scripts_dir, "mustache_circos.py"),
		run_circos = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/run_circos.py",
		outdir = absjoin("07.loops", "02.mustache", "{sample}", "{sample}_circos"),
		outpfix = absjoin("07.loops", "02.mustache", "{sample}", "{sample}_circos", "{sample}"),
		resolution = sorted(map(int, get_reslu("loop")))[0],
	shell:
		"export PATH=/work/frasergen/PUB/software/perl/perl-5.34.0/bin:$PATH ;\n"
		"{params.python3} {params.mustache_circos} --chromsize {input.reflen}"
		"    --resolution {params.resolution} --gene_annot {input.genes}"
		"    --input_loops {input.loop_combine} --outpfix {params.outpfix} ;\n"	

		"{params.python3} {params.run_circos} --method cis --karyotype {output.karyotype}"
		"    --gene_density {output.gene_density} --link_density {output.link_density}"
		"    --color_links {output.color_links} --outdir {params.outdir} ;\n"

def mustache_files(splname):
	file_list = list()
	for reslu in get_reslu("loop"):
		mustache = absjoin("07.loops", "02.mustache", splname, f"{splname}_{reslu}.mustache")
		file_list.append(mustache)
	return file_list
		

rule mustache_combine:
	input:
		genes = absjoin("01.ref", "genome_gene.bed"),
		loop_files = lambda wildcards: mustache_files(wildcards.sample),
	output:
		loop_combine = absjoin("07.loops", "02.mustache", "{sample}", "{sample}.combine.mustache"),
		loop_bedpe = absjoin("07.loops", "02.mustache", "{sample}", "{sample}.combine.bedpe"),
	params:
		python3 = software.get("python3"),
		loop_merge = absjoin(scripts_dir, "loop_merge.py"),
		#loop_bridge = absjoin(scripts_dir, "loop_bridge.py"),
		loop_bridge = "/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/loop_bridge.py",
		outpfix = absjoin("07.loops", "02.mustache", "{sample}", "{sample}_annot"),
	shell:
		"{params.python3} {params.loop_merge} --loop_files {input.loop_files} --outfile {output.loop_combine} ;\n"	
		
		"""
		awk -v OFS="\\t" 'NR>1{{print $1,$2,$3,$4,$5,$6}}' {output.loop_combine} > {output.loop_bedpe} ;\n
		"""

		"{params.python3} {params.loop_bridge} --loops_bedpe {output.loop_bedpe}"
		"     --gene_annot {input.genes} --annot_type bed --outpfix {params.outpfix} ;\n"

rule mustache_loop:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		loop_cool = absjoin("07.loops", "01.prepare", "{sample}_{loop_reslu}.cool"),
	output:
		loop_table = absjoin("07.loops", "02.mustache", "{sample}", "{sample}_{loop_reslu}.mustache"), 
	params:
		activate = software.get("activate"), 
	threads: 8
	shell:
		"source {params.activate} mustache ;\n"
		
		"mustache --pThreshold 0.1 --sparsityThreshold 0.88 --octaves 2 "
		"    --file {input.loop_cool} --outfile {output.loop_table}"
		"    --resolution {wildcards.loop_reslu} --chromosomeSize {input.reflen}"
		"    --processes {threads} ;\n"

rule loop_prepare:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		sample_pair = absjoin("03.align", "{sample}.merged.pair.gz"),
	output:
		loop_cool = absjoin("07.loops", "01.prepare", "{sample}_{loop_reslu}.cool"),
	params:
		#cooler = software.get("cooler"),
		cooler = "/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/cooler",
	shell:
		"{params.cooler} cload pairs {input.reflen}:{wildcards.loop_reslu}"
		"    {input.sample_pair} {output.loop_cool} -c1 2 -p1 3 -c2 4 -p2 5 ;\n"
		
		"{params.cooler} balance {output.loop_cool} ;\n"

rule cworld_insulation:
	input:
		inref = absjoin("01.ref", "genome.fa"),
		reflen = absjoin("01.ref", "genome.chrsize"),
		genes = absjoin("01.ref", "genome_gene.bed"),
		motif_db = config.get("jaspar_db"),
		tad_cool = absjoin("06.tad", "01.prepare", "{sample}_{tad_reslu}.tad.cool"),
		tad_basic_cool = absjoin("06.tad", "01.prepare", "{sample}.basic.cool"),
	output:
		interior_boundray = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}", 
							"boundary_interior.bed"),
		gc_gene = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}", 
							"boundary_interior.gc.gene-density"),
		is_bedgraph = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}.insulation.bedgraph"),
		is_bw = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}.score.bw"),
		fill_bw = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}",
							"{sample}_{tad_reslu}_fillnan.bw"),
		domain_bed = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}", 
							"TAD_genome.bed"),
		domain_ata = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}.ata"),
		top5_motif_png = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}",
								"{sample}.{tad_reslu}_boundary_motif5.png"),
		top5_motif_pdf = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}",
								"{sample}.{tad_reslu}_boundary_motif5.pdf"),

	params:
		genome_motif = config.get("genome_motif"),
		perl = software.get("perl"),
		python = software.get("python3"),
		cworld = software.get("cworld"),
		coolpup = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/cool2ata.py",
		rscript = software.get("Rscript"),
		ceqlogo = software.get("ceqlogo"),
		bedtools = software.get("bedtools"),
		hicexplorer = software.get("hicexplorer"),
		bdg2bw = absjoin(src_dir, "bedGraphToBigWig"),
		hic3_util = absjoin(scripts_dir, "hic3_util.py"),
		figure_util = absjoin(scripts_dir, "figure_util.py"),
		tad_display = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/tad_display.py",
		gc_gene_plot = absjoin(scripts_dir, "tad_gc_gene_plot.r"),
		outdir = absjoin("06.tad", "04.cworld", "{sample}.{tad_reslu}"),
		cooler_insulation = absjoin(scripts_dir, "cooler_insulation.py"),
		tad_boundary_motif = absjoin(scripts_dir, "tad_boundary_motif.py"),
	threads: 8
	shell:
		"{params.python} {params.cooler_insulation} --incool {input.tad_cool}"
		"    --cworld {params.cworld} --perl {params.perl} --outdir {params.outdir} ;\n"
		
		"{params.python} {params.hic3_util} bins_gcgene --inref {input.inref}"
		"    --gene-bed {input.genes} --region-bed {output.interior_boundray}"
		"    --outfile {output.gc_gene} ;\n"

		"{params.rscript} {params.gc_gene_plot} {output.gc_gene} {params.outdir} ;\n"

		""" for i in $(cut -f 1 {input.reflen}); do \
				cat {params.outdir}/${{i}}/*insulation.bedGraph; \
			done | awk '$4!="NA"' | {params.bedtools} sort > {output.is_bedgraph}; \n\
		"""
		
		"{params.bdg2bw} {output.is_bedgraph} {input.reflen} {output.is_bw} ;\n"
		"res=$[2*{wildcards.tad_reslu}] && bigwigCompare -b1 {output.is_bw} -b2 {output.is_bw} -bs $res -o {output.fill_bw} -p 10 --operation mean ;\n"
		"{params.python} {params.tad_display} --inref {input.inref}"
		"    --binsize {wildcards.tad_reslu} --cool-matrix {input.tad_cool}"
		"    --is-bw {output.fill_bw} --tad-bed {output.domain_bed}"
		"    --outdir {params.outdir}/TAD_Display ;\n"
		
		"{params.coolpup} tad {input.tad_basic_cool} {output.domain_bed} {output.domain_ata} ;\n"	

		# Do not use --minsize/--maxsize to filt TAD ! because sometime it comes out an error 
		#"{params.coolpup} --local --rescale --minsize 40000 --maxsize 3000000 --n_proc {threads}"
		#"     --outname {output.domain_ata} --seed 100 {input.tad_basic_cool} {output.domain_bed} ;\n"

		"{params.python} {params.tad_boundary_motif} --boundary_interior {output.interior_boundray}"
		"     --genome_motif {params.genome_motif} --ceqlogo {params.ceqlogo} --outdir {params.outdir}"
		"     --motif_db {input.motif_db} ;\n"

		"{params.python} {params.figure_util} combine_figure --out-image {output.top5_motif_png}"
		"     --ncol 1 {params.outdir}/Motif1.temp.png {params.outdir}/Motif2.temp.png"
		"     {params.outdir}/Motif3.temp.png {params.outdir}/Motif4.temp.png {params.outdir}/Motif5.temp.png ;\n"

		"{params.python} {params.figure_util} combine_figure --out-image {output.top5_motif_pdf}"
		"     --ncol 1 {params.outdir}/Motif1.temp.png {params.outdir}/Motif2.temp.png"
		"     {params.outdir}/Motif3.temp.png {params.outdir}/Motif4.temp.png {params.outdir}/Motif5.temp.png ;\n"

		"rm {params.outdir}/Motif*.temp.png ;\n"


rule run_hicfindtad:
	input:
		inref = absjoin("01.ref", "genome.fa"),
		reflen = absjoin("01.ref", "genome.chrsize"),
		genes = absjoin("01.ref", "genome_gene.bed"),
		motif_db = config.get("jaspar_db"),
		tad_cool = absjoin("06.tad", "01.prepare", "{sample}_{tad_reslu}.tad.cool"),
		tad_basic_cool = absjoin("06.tad", "01.prepare", "{sample}.basic.cool"),
	output:
		domain_info = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}", 
						"{sample}.{tad_reslu}_domains.bed"),
		domain_bed = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}", 
							"TAD_genome.bed"),
		boundary_info = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}_boundaries.bed"),
		interior_boundray = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}", 
							"boundary_interior.bed"),
		gc_gene = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}", 
							"boundary_interior.gc.gene-density"),
		domain_bedg = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}_score.bedgraph"),
		domain_bedg_sort = temp(absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}_score.bedgraph.sort")),	
		domain_bw = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}.score.bw"),	
		fill_bw = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}_fillnan.bw"),	
		domain_ata = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}",
							"{sample}.{tad_reslu}.ata"),
		top5_motif_png = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}",
								"{sample}.{tad_reslu}_boundary_motif5.png"),
		top5_motif_pdf = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}",
								"{sample}.{tad_reslu}_boundary_motif5.pdf"),
	params:
		genome_motif = config.get("genome_motif"),
		python = software.get("python3"),
		rscript = software.get("Rscript"),
		coolpup = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/cool2ata.py",
		ceqlogo = software.get("ceqlogo"),
		bedtools = software.get("bedtools"),
		hicexplorer = software.get("hicexplorer"),
		bdg2bw = absjoin(src_dir, "bedGraphToBigWig"),
		hic3_util = absjoin(scripts_dir, "hic3_util.py"),
		figure_util = absjoin(scripts_dir, "figure_util.py"),
		tad_display = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/tad_display.py",
		gc_gene_plot = absjoin(scripts_dir, "tad_gc_gene_plot.r"),
		post_hicfindtad = absjoin(scripts_dir, "post_hicfindtad.py"),
		tad_boundary_motif = absjoin(scripts_dir, "tad_boundary_motif.py"),
		outdir = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}"),
		hicFindTADs = absjoin(dirname(software.get("hicexplorer")), "hicFindTADs"),
		outpfix = absjoin("06.tad", "03.hicFindTAD", "{sample}.{tad_reslu}", "{sample}.{tad_reslu}"),
	threads: 8
	shell:
		"{params.hicFindTADs} --matrix {input.tad_cool} --outPrefix {params.outpfix}"
		"    --numberOfProcessors {threads} --correctForMultipleTesting fdr ;\n"
		
		"cut -f1-3 {output.domain_info} > {output.domain_bed} ;\n"

		"{params.python} {params.post_hicfindtad} --chromsize {input.reflen}"
		"    --tad-bed {output.domain_bed} --binsize {wildcards.tad_reslu}"
		"    --boundary {output.boundary_info} --outdir {params.outdir} ;\n"
	
		"{params.python} {params.hic3_util} bins_gcgene --inref {input.inref}" 
		"    --gene-bed {input.genes} --region-bed {output.interior_boundray}"
		"    --outfile {output.gc_gene} ;\n"

		"{params.rscript} {params.gc_gene_plot} {output.gc_gene} {params.outdir} ;\n"

		"{params.bedtools} sort -i {output.domain_bedg} > {output.domain_bedg_sort} ;\n"
		
		"{params.bdg2bw} {output.domain_bedg_sort} {input.reflen} {output.domain_bw} ;\n"
		"res=$[2*{wildcards.tad_reslu}] && bigwigCompare -b1 {output.domain_bw} -b2 {output.domain_bw} -bs $res -o {output.fill_bw} -p 10 --operation mean ;\n"
		"{params.python} {params.tad_display} --inref {input.inref}"
		"    --binsize {wildcards.tad_reslu} --cool-matrix {input.tad_cool}"
		"    --is-bw {output.fill_bw} --tad-bed {output.domain_info}"
		"    --outdir {params.outdir}/TAD_Display ;\n"

		"{params.coolpup} tad {input.tad_basic_cool} {output.domain_bed} {output.domain_ata} ;\n"

		"{params.python} {params.tad_boundary_motif} --boundary_interior {output.interior_boundray}"
		"     --genome_motif {params.genome_motif} --ceqlogo {params.ceqlogo} --outdir {params.outdir}"
		"     --motif_db {input.motif_db} ;\n"

		"{params.python} {params.figure_util} combine_figure --out-image {output.top5_motif_png}"
		"     --ncol 1 {params.outdir}/Motif1.temp.png {params.outdir}/Motif2.temp.png"
		"     {params.outdir}/Motif3.temp.png {params.outdir}/Motif4.temp.png {params.outdir}/Motif5.temp.png ;\n"

		"{params.python} {params.figure_util} combine_figure --out-image {output.top5_motif_pdf}"
		"     --ncol 1 {params.outdir}/Motif1.temp.png {params.outdir}/Motif2.temp.png"
		"     {params.outdir}/Motif3.temp.png {params.outdir}/Motif4.temp.png {params.outdir}/Motif5.temp.png ;\n"

		"rm {params.outdir}/Motif*.temp.png ;\n"

		
rule diamond_insulation:
	input:
		inref = absjoin("01.ref", "genome.fa"),
		genes = absjoin("01.ref", "genome_gene.bed"),
		reflen = absjoin("01.ref", "genome.chrsize"),
		motif_db = config.get("jaspar_db"),
		tad_cool = absjoin("06.tad", "01.prepare", "{sample}_{tad_reslu}.tad.cool"),
		tad_basic_cool = absjoin("06.tad", "01.prepare", "{sample}.basic.cool"),
	output:
		diamond = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}", "{sample}.{tad_reslu}.tsv"),
		diamond_domain = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}", "TAD_genome.bed"),
		diamond_ib = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}", "boundary_interior.bed"),
		diamond_gcgene = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}", 
								"boundary_interior.gc.gene-density"),
		diamond_bw = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}", 
								"{sample}.{tad_reslu}.score.bw"),
		fill_bw = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}", 
								"{sample}.{tad_reslu}_fillnan.bw"),
		diamond_ata = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}",
								"{sample}.{tad_reslu}.ata"),
		top5_motif_png = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}",
								"{sample}.{tad_reslu}_boundary_motif5.png"),
		top5_motif_pdf = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}",
								"{sample}.{tad_reslu}_boundary_motif5.pdf"),
	params:
		genome_motif = config.get("genome_motif"),
		python = software.get("python3"),
		rscript = software.get("Rscript"),
		coolpup = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/cool2ata.py",
		ceqlogo = software.get("ceqlogo"),
		cooltools = software.get("cooltools"),
		hicexplorer = software.get("hicexplorer"),
		hic3_util = absjoin(scripts_dir, "hic3_util.py"),
		tad_display = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/tad_display.py",
		figure_util = absjoin(scripts_dir, "figure_util.py"),
		gc_gene_plot = absjoin(scripts_dir, "tad_gc_gene_plot.r"),
		outdir = absjoin("06.tad", "02.cooltools", "{sample}.{tad_reslu}"),
		post_diamond = absjoin(scripts_dir, "post_diamond.py"),
		tad_boundary_motif = absjoin(scripts_dir, "tad_boundary_motif.py"),
		window_size = lambda wildcards: int(wildcards.tad_reslu) * 10,
	threads: 8
	shell:
		"{params.cooltools} insulation --ignore-diags 2 --threshold 0.1"
		"    --output {output.diamond} --min-dist-bad-bin 5 --append-raw-scores"
		"    --bigwig {input.tad_cool} {params.window_size} ;\n"

		"mv {params.outdir}/{wildcards.sample}.{wildcards.tad_reslu}.tsv.*.bw"
		"    {output.diamond_bw} ;\n"

		"{params.python} {params.post_diamond} --chrsize {input.reflen}" 
		"    --diamond_insulation {output.diamond}"
		"    --sample_name {wildcards.sample}.{wildcards.tad_reslu}"
		"    --resolution {wildcards.tad_reslu} --outdir {params.outdir} ;\n"
		
		"{params.python} {params.hic3_util} bins_gcgene --inref {input.inref}"
		"    --gene-bed {input.genes} --region-bed {output.diamond_ib}"
		"    --outfile {output.diamond_gcgene} ;\n"

		"{params.rscript} {params.gc_gene_plot} {output.diamond_gcgene} {params.outdir} ;\n"

		"res=$[2*{wildcards.tad_reslu}] && bigwigCompare -b1 {output.diamond_bw} -b2 {output.diamond_bw} -bs $res -o {output.fill_bw} -p 10 --operation mean ;\n"
		"{params.python} {params.tad_display} --inref {input.inref}"
		"     --binsize {wildcards.tad_reslu} --cool-matrix {input.tad_cool}"
		"     --is-bw {output.fill_bw} --tad-bed {output.diamond_domain}"
		"     --outdir {params.outdir}/TAD_Display ;\n"

		"{params.coolpup} tad {input.tad_basic_cool} {output.diamond_domain} {output.diamond_ata} ;\n"

		"{params.python} {params.tad_boundary_motif} --boundary_interior {output.diamond_ib}"
		"     --genome_motif {params.genome_motif} --ceqlogo {params.ceqlogo}"
		"     --motif_db {input.motif_db} --outdir {params.outdir} ;\n"

		"{params.python} {params.figure_util} combine_figure --out-image {output.top5_motif_png}"
		"     --ncol 1 {params.outdir}/Motif1.temp.png {params.outdir}/Motif2.temp.png"
		"     {params.outdir}/Motif3.temp.png {params.outdir}/Motif4.temp.png {params.outdir}/Motif5.temp.png ;\n"

		"{params.python} {params.figure_util} combine_figure --out-image {output.top5_motif_pdf}"
		"     --ncol 1 {params.outdir}/Motif1.temp.png {params.outdir}/Motif2.temp.png"
		"     {params.outdir}/Motif3.temp.png {params.outdir}/Motif4.temp.png {params.outdir}/Motif5.temp.png ;\n"
	
		"rm {params.outdir}/Motif*.temp.png ;\n"

rule tad_prepare:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		sample_pair = absjoin("03.align", "{sample}.merged.pair.gz"),
	output:
		tad_cool = absjoin("06.tad", "01.prepare", "{sample}_{tad_reslu}.tad.cool"),
	params:
		#cooler = software.get("cooler"),
		cooler = "/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/cooler",
	shell:
		"{params.cooler} cload pairs {input.reflen}:{wildcards.tad_reslu}"
		"    {input.sample_pair} {output.tad_cool} -c1 2 -p1 3 -c2 4 -p2 5 ;\n"
		
		"{params.cooler} balance {output.tad_cool} ;\n"

rule tad_prepare_basic:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		sample_pair = absjoin("03.align", "{sample}.merged.pair.gz"),
	output:
		tad_basic_cool = absjoin("06.tad", "01.prepare", "{sample}.basic.cool"),
	params:
		basic_reslu = config.get("resolution").get("basic"),
		#cooler = software.get("cooler"),
		cooler = "/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/cooler",
		logfile = absjoin("04.matrix", "work.sh"),
	threads: 8 
	shell:
		"{params.cooler} cload pairs {input.reflen}:{params.basic_reslu}"
		"    {input.sample_pair} {output.tad_basic_cool} -c1 2 -p1 3 -c2 4 -p2 5 ;\n"
		
		"{params.cooler} balance --nproc {threads} --max-iters 500 {output.tad_basic_cool} ;\n"

rule run_cscore:
	input:
		inref = absjoin("01.ref", "genome.fa"),
		reflen = absjoin("01.ref", "genome.chrsize"),
		genes = absjoin("01.ref", "genome_gene.bed"),
		contact_file = absjoin("03.align", "{sample}.contact"),
		cptmt_cool = absjoin("05.compartment", "01.prepare", "{sample}_{cptmt_reslu}.cool"),
		cscore_bin = absjoin("05.compartment", "01.prepare", "genome_cscore_bin_{cptmt_reslu}.bed"),
	output:
		cscore_cptmt = absjoin("05.compartment", "03.Cscore", "{sample}_{cptmt_reslu}", 
							"{sample}_{cptmt_reslu}_cscore.bedgraph"),
		cscore_cptmt_bw = absjoin("05.compartment", "03.Cscore", "{sample}_{cptmt_reslu}", 
							"{sample}_{cptmt_reslu}.bw"),
		cscore_cptmt_gcgene = absjoin("05.compartment", "03.Cscore", "{sample}_{cptmt_reslu}", 
							"{sample}_{cptmt_reslu}.gc_gene.bed"),
	params:
		splres = "{sample}_{cptmt_reslu}",
		cscore = software.get("cscore"),
		perl = software.get("perl"),
		python = software.get("python3"),
		rscript = software.get("Rscript"),
		cptmt_merge = absjoin(scripts_dir, "compartment_merge.pl"),
		cptmt_gcplot = absjoin(scripts_dir, "cptmt_gc_gene_plot.r"),
		bdg2bw = absjoin(src_dir, "bedGraphToBigWig"),
		hic3_util = absjoin(scripts_dir, "hic3_util.py"),
		outdir = absjoin("05.compartment", "03.Cscore", "{sample}_{cptmt_reslu}"),
		outpfix = absjoin("05.compartment", "03.Cscore", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}"),
		cscore_cptmt_hm_prefix = absjoin("05.compartment", "03.Cscore", "{sample}_{cptmt_reslu}", 
							"{sample}_{cptmt_reslu}"),
		logfile = absjoin("05.compartment", "03.Cscore", "work.sh"),
		hicPlotMatrix = absjoin(dirname(software.get("hicexplorer")), "hicPlotMatrix"),
	shell:
		"{params.cscore} {input.cscore_bin} {input.contact_file} {params.outpfix} 8 {wildcards.cptmt_reslu} ;\n"
		"echo {params.cscore} {input.cscore_bin} {input.contact_file} {params.outpfix}"
		"    8 {wildcards.cptmt_reslu} >> {params.logfile} ;\n\n"

		"{params.bdg2bw} {output.cscore_cptmt} {input.reflen} {output.cscore_cptmt_bw} ;\n"
		"echo {params.bdg2bw} {output.cscore_cptmt} {input.reflen} {output.cscore_cptmt_bw} >> {params.logfile} ;\n\n"

		"for i in $(cut -f 1 {input.reflen}); do "
		"    {params.hicPlotMatrix} --matrix {input.cptmt_cool}"
		"        --outFileName {params.cscore_cptmt_hm_prefix}.${{i}}.bychr.png"
		"        --fontsize 14 --log --dpi 800 --region ${{i}} --bigwig {output.cscore_cptmt_bw}"
		"        --vMinBigwig -1 --vMaxBigwig 1; done ;\n"

		"{params.python} {params.hic3_util} bins_gcgene --inref {input.inref} --gene-bed {input.genes}"
		"        --region-bed {output.cscore_cptmt} --outfile {output.cscore_cptmt_gcgene} ;\n"

		"{params.rscript} {params.cptmt_gcplot} {output.cscore_cptmt_gcgene} {params.outdir} ;\n"
	
		"{params.perl} {params.cptmt_merge} -i {output.cscore_cptmt_gcgene} -o {params.outdir} -s {params.splres} ;\n"

rule csocre_contact:
	input:
		sample_pair = absjoin("03.align", "{sample}.merged.pair.gz"),
	output:
		contact_file = temp(absjoin("03.align", "{sample}.contact")),
	params:
		bgzip = software.get("bgzip"),
	threads: 12
	shell:
		"""{params.bgzip} -dc -@ {threads} {input.sample_pair} |\
			  awk -v OFS="\\t" '$0!~/^#/{{print $1,$2,$3,$6,$4,$5,$7}}' > {output.contact_file} ;\n"""

rule run_pca:
	input:
		inref = absjoin("01.ref", "genome.fa"),
		reflen = absjoin("01.ref", "genome.chrsize"),
		genes = absjoin("01.ref", "genome_gene.bed"),
		genome_gene = absjoin("05.compartment", "01.prepare", "genome_genecov_{cptmt_reslu}.bed"),
		cptmt_cool = absjoin("05.compartment", "01.prepare", "{sample}_{cptmt_reslu}.cool"),
	output:
		exp_table = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}.exp.tsv"),
		all_eigs = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}.eigs.cis.vecs.tsv"),
		cptmt_bdg = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}.bedgraph"),
		cptmt_output = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}.compartment.xls"),
		cptmt_bw = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}.bw"),
		cptmt_gcgene = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}.gc_gene.bed"),
		hicpca_out = temp(absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}.hicpca.bed")),
		pearson_hm = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}.phm"),
	params:
		splres = "{sample}_{cptmt_reslu}",
		perl = software.get("perl"),
		python = software.get("python3"),
		bedtools = software.get("bedtools"),
		rscript = software.get("Rscript"),
		cooltools = software.get("cooltools"),
		hic3_util = absjoin(scripts_dir, "hic3_util.py"),
		bedGraphToBigWig = software.get("bedGraphToBigWig"),
		figure_util = absjoin(scripts_dir, "figure_util.py"),
		compartment_plot = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/compartment_plot.py",
		# old pearson plot using plotgardener is deprecated
		# pearson_plot = absjoin(scripts_dir, "pearson_plot.r"),
		#rscript_plotgardener = software.get("Rscript_plotgardener"),
		cptmt_merge = absjoin(scripts_dir, "compartment_merge.pl"),
		cptmt_gcplot = absjoin(scripts_dir, "cptmt_gc_gene_plot.r"),
		hicPCA = absjoin(dirname(software.get("hicexplorer")), "hicPCA"),
		outdir = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}"),
		outpfix = absjoin("05.compartment", "02.PCA", "{sample}_{cptmt_reslu}", "{sample}_{cptmt_reslu}"),
	shell:
		"{params.cooltools} expected-cis --nproc 32 --output {output.exp_table} {input.cptmt_cool} ;\n"
		
		"{params.cooltools} eigs-cis --phasing-track {input.genome_gene} --out-prefix {params.outpfix}.eigs"
		"    {input.cptmt_cool} ;\n"
		
		"{params.cooltools} saddle --contact-type cis --strength --out-prefix {params.outpfix}.saddle"
		"    --qrange 0.02 0.98 --fig pdf {input.cptmt_cool} {params.outpfix}.eigs.cis.vecs.tsv {output.exp_table} ;\n"
		
		"{params.cooltools} saddle --contact-type cis --strength --out-prefix {params.outpfix}.saddle"
		"    --qrange 0.02 0.98 --fig png {input.cptmt_cool} {params.outpfix}.eigs.cis.vecs.tsv {output.exp_table} ;\n"

		"""awk -v OFS="\\t" 'NR>1{{if (NF==7) {{print $1,$2,$3,$5}} else {{print $1,$2,$3,0}}}}' \
			{params.outpfix}.eigs.cis.vecs.tsv | {params.bedtools} sort > {output.cptmt_bdg} ;\n """

		"{params.python} {params.hic3_util} bins_gcgene --inref {input.inref} --gene-bed {input.genes}"
		"    --region-bed {output.cptmt_bdg} --outfile {output.cptmt_gcgene} ;\n"

		"{params.rscript} {params.cptmt_gcplot} {output.cptmt_gcgene} {params.outdir} ;\n"

		"{params.perl} {params.cptmt_merge} -i {output.cptmt_gcgene} -o {params.outdir} -s {params.splres} ;\n"
		
		"{params.hicPCA} --matrix {input.cptmt_cool} --outputFileName {output.hicpca_out} --whichEigenvectors 1"
		"    --pearsonMatrix {output.pearson_hm} ;\n"
		
		"{params.bedGraphToBigWig} {output.cptmt_bdg} {input.reflen} {output.cptmt_bw} ;\n"
		
		
		"{params.python} {params.compartment_plot} --chromsize {input.reflen} --control_name {params.splres}"
		"    --control_cool {input.cptmt_cool} --control_cptmt_gcgene {output.cptmt_gcgene} --pc1_column 4"
		"    --outdir {params.outpfix} ;\n"

		# old pearson plot using plotgardener is deprecated
		#"{params.rscript_plotgardener} {params.pearson_plot} --pearson_hm {output.pearson_hm}"
		#"    --cptmt_bw {output.cptmt_bw} --outdir {params.outdir} --splname {params.splres} ;\n"

		#"for i in $(cut -f 1 {input.reflen}); do"
		#"    {params.python} {params.figure_util} pdf2png --inpdf {params.outpfix}.${{i}}.bychr.pdf "
		#"    --outpng {params.outpfix}.${{i}}.bychr.png; done ;\n"

		"""awk 'BEGIN{{print "chrom\\tstart\\tend\\tPC1\\tCompartment"}} \
			{{if ($4>0){{print $0"\\tA"}} else if ($4<0) {{print $0"\\tB"}} else {{print $0"\\t-"}}}}' \
			{output.cptmt_bdg} > {output.cptmt_output} ;\n """

rule cptmt_prepare_cool:
	input:
		basic_cool = absjoin("04.matrix", "{sample}.basic.cool"),
	output:
		cptmt_cool = absjoin("05.compartment", "01.prepare", "{sample}_{cptmt_reslu}.cool"),
	params:
		#cooler = software.get("cooler"),
		cooler = "/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/cooler",
		basic_reslu = config.get("resolution").get("basic"),
		logfile = absjoin("05.compartment", "01.prepare", "work.sh"),
	threads: 8
	shell:
		"if [ {wildcards.cptmt_reslu} -eq {params.basic_reslu} ]; then"
		"    cp {input.basic_cool} {output.cptmt_cool}; else "
		"{params.cooler} coarsen --nproc {threads} --factor $[{wildcards.cptmt_reslu}/{params.basic_reslu}]"
		"    --out {output.cptmt_cool} {input.basic_cool} ; fi \n"
		"echo {params.cooler} coarsen --nproc {threads} --factor $[{wildcards.cptmt_reslu}/{params.basic_reslu}]"
		"    --out {output.cptmt_cool} {input.basic_cool} >> {params.logfile} ;\n\n"
		
		"{params.cooler} balance {output.cptmt_cool} ;\n"
		"echo {params.cooler} balance {output.cptmt_cool} >> {params.logfile} ;\n\n"

rule cptmt_prepare_bin:
	input:
		inref = absjoin("01.ref", "genome.fa"),
		genes = absjoin("01.ref", "genome_gene.bed"),
		reflen = absjoin("01.ref", "genome.chrsize"),
	output:
		genome_bin = absjoin("05.compartment", "01.prepare", "genome_bin_{cptmt_reslu}.bed"),
		genome_gc = absjoin("05.compartment", "01.prepare", "genome_gc_{cptmt_reslu}.bed"),
		genome_gene = absjoin("05.compartment", "01.prepare", "genome_genecov_{cptmt_reslu}.bed"),
		cscore_bin = absjoin("05.compartment", "01.prepare", "genome_cscore_bin_{cptmt_reslu}.bed"),
	params:
		python = software.get("python3"),
		hic3_util = absjoin(scripts_dir, "hic3_util.py"),
		cooltools = software.get("cooltools"),
		equal_len_bed = absjoin(dirname(software.get("cscore")), "generateEqualLengthBed"),		
	shell:
		"{params.cooltools} genome binnify {input.reflen} {wildcards.cptmt_reslu} | grep chr > {output.genome_bin} ;\n"
	
		"{params.cooltools} genome gc {output.genome_bin} {input.inref} > {output.genome_gc} ;\n"
	
		"{params.equal_len_bed} {input.reflen} {output.cscore_bin} {wildcards.cptmt_reslu} ;\n"

		"{params.python} {params.hic3_util} bin_genecov --bin-bed {output.genome_bin}"
		"    --gene-bed {input.genes} --outfile {output.genome_gene} ;\n"

rule cool_matrix:
	input:
		absjoin("04.matrix", "{sample}.{reslu}.cool"),
	output:
		absjoin("04.matrix", "{sample}_{reslu}", "{sample}.zScore.matrix.gz"),
	run:
		parse_genome_cool(*input,output[0].replace("zScore.",""))
rule cool_matrix_bychr:
	input:
		absjoin("05.compartment", "01.prepare", "{sample}_{cptmt_reslu}.cool"),
	output:
		absjoin("04.matrix","{sample}_{cptmt_reslu}","bychr","{sample}.chr1.zScore.matrix.gz")
	run:
		parse_chrom_cool(*input,output[0].replace("chr1.zScore.",""))

rule heatmap_distance:
	input:
		basic_cool = absjoin("04.matrix", "{sample}.basic.cool"),
	output:
		genome_cool = absjoin("04.matrix", "{sample}.{genome_reslu}.cool"),
		genome_heatmap = absjoin("04.matrix", "{sample}_{genome_reslu}.genome.png"),
		chr1_heatmap = absjoin("04.matrix", "{sample}_{genome_reslu}", "{sample}_{genome_reslu}.chr1.png"),
		ps_curve = absjoin("04.matrix", "{sample}.ps_curve.{genome_reslu}.png"),
		ps_curve_data = absjoin("04.matrix", "{sample}.ps_curve.{genome_reslu}.txt"),
	params:
		python3 = software.get("python3"),
		cool_heatmap = absjoin(scripts_dir, "cool_heatmap.py"),
		basic_reslu = config.get("resolution").get("basic"),
		#cooler = software.get("cooler"),
		cooler = "/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/cooler",
		logfile = absjoin("04.matrix", "work.sh"),
		genome_outpfix = absjoin("04.matrix", "{sample}_{genome_reslu}"),
		chrom_outpfix = absjoin("04.matrix", "{sample}_{genome_reslu}", "{sample}_{genome_reslu}"),
		hicPlotDistVsCounts = absjoin(dirname(software.get("hicexplorer")), "hicPlotDistVsCounts"),
	threads: 8
	shell:
		"if [ {wildcards.genome_reslu} -eq {params.basic_reslu} ]; then"
		"    cp {input.basic_cool} {output.genome_cool}; else "
		"{params.cooler} coarsen --nproc {threads} --factor $[{wildcards.genome_reslu}/{params.basic_reslu}]"
		"    --out {output.genome_cool} {input.basic_cool} ; fi \n"
	
		"{params.cooler} balance {output.genome_cool} ;\n"
		"echo {params.cooler} balance {output.genome_cool} >> {params.logfile} ;\n\n"

		"{params.python3} {params.cool_heatmap} --incool {output.genome_cool} --region genome"
		"    --outpfix {params.genome_outpfix} --vmax 0.001 ;\n"
		"echo {params.python3} {params.cool_heatmap} --incool {output.genome_cool} --region genome"
		"    --outpfix {params.genome_outpfix} --vmax 0.001 >> {params.logfile} ;\n\n"
		
		"{params.python3} {params.cool_heatmap} --incool {output.genome_cool} --region chrom"
		"     --outpfix {params.chrom_outpfix} --vmax 0.001 ;\n"
		"echo {params.python3} {params.cool_heatmap} --incool {output.genome_cool} --region chrom"
		"     --outpfix {params.chrom_outpfix} --vmax 0.001 >> {params.logfile} ;\n\n"
	
		"{params.hicPlotDistVsCounts} --matrices {output.genome_cool} --labels {wildcards.sample}"
		"    --plotFile {output.ps_curve} --outFileData {output.ps_curve_data} --plotsize 10 7 ;\n"
		"echo {params.hicPlotDistVsCounts} --matrices {output.genome_cool} --labels {wildcards.sample}"
		"    --plotFile {output.ps_curve} --outFileData {output.ps_curve_data} --plotsize10 7 >> {params.logfile} ;\n\n"

rule basic_cool:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		sample_pair = absjoin("03.align", "{sample}.merged.pair.gz"),
	output:
		basic_cool = absjoin("04.matrix", "{sample}.basic.cool"),
	params:
		basic_reslu = config.get("resolution").get("basic"),
		#cooler = software.get("cooler"),
		cooler = "/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/cooler",
		logfile = absjoin("04.matrix", "work.sh"),
	shell:
		"{params.cooler} cload pairs {input.reflen}:{params.basic_reslu}"
		"    {input.sample_pair} {output.basic_cool} -c1 2 -p1 3 -c2 4 -p2 5 ;\n"
		"echo {params.cooler} cload pairs {input.reflen}:{params.basic_reslu}"
		"    {input.sample_pair} {output.basic_cool} -c1 2 -p1 3 -c2 4 -p2 5 > {params.logfile} ;\n\n"

rule pair_stat:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		sample_pair = absjoin("03.align", "{sample}.merged.pair.gz"),
	output:
		pairqc_stat = absjoin("03.align", "{sample}_report", "{sample}.cis_to_trans.out"),
		pairqc_fig = absjoin("03.align", "{sample}_report", "plots", "{sample}.proportion.png"),
	params:
		python = software.get("python3"),
		pairix = software.get("pairix"),
		pairqc = software.get("pairqc"),
		rscript = software.get("Rscript"),
		outpfix = absjoin("03.align", "{sample}"),
		pairqc_dir = absjoin("03.align", "{sample}_report"),
		pairqc_plot = absjoin(dirname(software.get("pairqc")), "plot.r"),
		max_dist = get_maxlog_dist(config.get("reference")),
	shell:
		"{params.pairix} -f -p pairs {input.sample_pair} ;\n"

		"{params.python} {params.pairqc} --pairs {input.sample_pair} --chrsize {input.reflen}"
		"    --input_type P --outdir_prefix {params.outpfix} --sample_name {wildcards.sample}"
		"    --max_logdistance {params.max_dist} ;\n"

		# 4 for 4-cut enzyme, eg. MboI
		"{params.rscript} {params.pairqc_plot} 4 {params.pairqc_dir} ;\n"

rule resolution_assess:
	input:
		infai = absjoin("01.ref", "genome.fa.fai"),
		sample_pair = absjoin("03.align", "{sample}.merged.pair.gz"),
	output:
		reslu_figure = absjoin("03.align", "{sample}.reslu.png"),
	params:
		python = software.get("python3"),
		outpfix = absjoin("03.align", "{sample}"),
		assess_reslu = absjoin(scripts_dir, "assess_reslu.py"),
	resources:
		mem_mb=300
	shell:
		"{params.python} {params.assess_reslu} --valid_pairs {input.sample_pair}"
		"    --infai {input.infai} --outpfix {params.outpfix} ;\n"

rule merge_pair:
	input:
		biorep_pair = lambda wildcards: [absjoin("03.align", f"{biorep}.nodup.pair.gz")
				for biorep in config[wildcards.sample].keys()],
	output:
		sample_pair = absjoin("03.align", "{sample}.merged.pair.gz"),
	params:
		bgzip = software.get("bgzip"),
		pairtools = software.get("pairtools"),
		logfile = absjoin("03.align", "work.sh"),
	threads: 10
	shell:
		"rm -rf {output.sample_pair} ;\n"

		"{params.pairtools} merge --output {output.sample_pair} --compress-program gzip {input.biorep_pair}"
		"    --cmd-out '{params.bgzip} -@ 12' --cmd-in '{params.bgzip} -dc -@ 12' ;\n"

		"echo {params.pairtools} merge --output {output.sample_pair} --compress-program gzip"
		"    {input.biorep_pair} --cmd-out \'{params.bgzip} -@ 10\' >> {params.logfile} ;\n\n"

rule compact_pair:
	input:
		nodup_pair = absjoin("03.align", "{dataset}.nodup.pair"), 
		dataset_contact = lambda wildcards: expand(absjoin("04.matrix", "reproduce_{reslu}", 
				wildcards.dataset, "refined.frag.{reslu}.contact.gz"), reslu=get_reslu("reproduce")),
	output:
		dspair_gz = absjoin("03.align", "{dataset}.nodup.pair.gz"),
	params:
		bgzip = software.get("bgzip"),
	shell:
		"echo {input.dataset_contact} ;\n"
		"{params.bgzip} -@ 12 {input.nodup_pair} ;\n"
	
rule combine_contact:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		contacts = reslu_contact("{reslu}"),
	output:
		logfile = absjoin("04.matrix", "reproduce_{reslu}", "work.sh"),
	params:
		datasets = " ".join(datasets),
		python = software.get("python3"),
		genomedisco = software.get("genomedisco"),
		outdir = absjoin("04.matrix", "reproduce_{reslu}"),
		reproduce_heatmap = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/reproduce_heatmap.py",
	run:
		script=rf"""
bedtools makewindows -g {input.reflen} -w {wildcards.reslu}|perl -alne 'print "$_\t".int(($F[1]+$F[2])/2)'|gzip > genome_bin200000.bed.gz
for i in {params.datasets}; do
echo -e $i"\t"{params.outdir}/$i/refined.frag.{wildcards.reslu}.contact.gz >> {params.outdir}/datasets_contact.fofn
for j in {params.datasets}; do
cat <<END >>script.sh
mkdir $i.vs.$j && cd $i.vs.$j
echo -e "$i\t{params.outdir}/$i/refined.frag.{wildcards.reslu}.contact.gz
$j\t{params.outdir}/$j/refined.frag.{wildcards.reslu}.contact.gz" > $i.vs.$j.contact
echo -e "$i\t$j" > $i.vs.$j.pair
{params.genomedisco} run_all --metadata_samples $i.vs.$j.contact --metadata_pairs $i.vs.$j.pair --bins ../genome_bin200000.bed.gz --outdir ./
cd ..
END
done
done

#{params.python} {params.reproduce_heatmap} --contact_fofn {params.outdir}/datasets_contact.fofn --chromsize {input.reflen} --genomedisco {params.genomedisco} --outdir {params.outdir} --resolution {wildcards.reslu} --labels {config["samplenames"]} --lens {{0}}
""".format(",".join(map(str,lens)))
		open(output.logfile,"a").write(script)

rule reproducibility:
	input:
		reflen = absjoin("01.ref", "genome.chrsize"),
		nodup_pair = absjoin("03.align", "{dataset}.nodup.pair"),
	output:
		dataset_cool = absjoin("04.matrix", "reproduce_{reslu}", "{dataset}.cool"),
		dataset_contact = absjoin("04.matrix", "reproduce_{reslu}", "{dataset}", 
						"refined.frag.{reslu}.contact.gz"),
	params:
		python = software.get("python3"),
		#cooler = software.get("cooler"),
		cooler = "/work/frasergen/backup/3d/project/Interactome/240590_human_20240717/03.HiC/cooler",
		outpfix = absjoin("04.matrix", "reproduce_{reslu}", "{dataset}"),
		cooler_fithic = absjoin(scripts_dir, "cooler_fithic2.py"),
	resources:
		mem_mb=200
	shell:
		"{params.cooler} cload pairs {input.reflen}:{wildcards.reslu} "
		"    {input.nodup_pair} {output.dataset_cool} -c1 2 -p1 3 -c2 4 -p2 5 ;\n"
		
		"{params.python} {params.cooler_fithic} -i {output.dataset_cool} -o {params.outpfix} ;\n"

rule aln_pair:
	input:
		inref = absjoin("01.ref", "genome.fa"),
		reflen = absjoin("01.ref", "genome.chrsize"),
		read1 = absjoin("02.reads", "{dataset}", "{dataset}_1P.fastq.gz"),
		read2 = absjoin("02.reads", "{dataset}", "{dataset}_2P.fastq.gz"),
	output:
		aln_bam = temp(absjoin("03.align", "{dataset}.sort.bam")), 
		raw_pair = temp(absjoin("03.align", "{dataset}.raw.pair")),
		sort_pair = temp(absjoin("03.align", "{dataset}.sort.pair")),
		nodup_pair = absjoin("03.align", "{dataset}.nodup.pair"),
		pair_stat = absjoin("03.align", "{dataset}.pair.stat"),
	params:
		bwa = software.get("bwa"),
		samtools = software.get("samtools"),
		pairtools = software.get("pairtools"),
		logfile = absjoin("03.align", "work.sh"),
	threads: 32
	shell:
		"{params.bwa} mem -SP5M -t {threads} {input.inref} {input.read1} {input.read2} |"
		"    {params.samtools} view --threads {threads} -Sb -> {output.aln_bam} ;\n"
		"echo {params.bwa} mem -SP5M -t {threads} {input.inref} {input.read1} {input.read2} \|"
		"    {params.samtools} view --threads {threads} -Sb -\> {output.aln_bam} > {params.logfile} ;\n\n"

		"{params.pairtools} parse -c {input.reflen} -o {output.raw_pair} --drop-sam {output.aln_bam} ;\n"
		"echo {params.pairtools} parse -c {input.reflen} -o {output.raw_pair} --drop-sam {output.aln_bam}"
		"    >> {params.logfile} ;\n\n"
	
		"{params.pairtools} sort --nproc {threads} -o {output.sort_pair} {output.raw_pair} ;\n"
		"echo {params.pairtools} sort --nproc {threads} -o {output.sort_pair} {output.raw_pair}"
		"    >> {params.logfile} ;\n\n"
	
		"{params.pairtools} dedup {output.sort_pair} --output {output.nodup_pair} ;\n"
		"echo {params.pairtools} dedup {output.sort_pair} --output {output.nodup_pair}"
		"    >> {params.logfile} ;\n\n"
		
		"{params.pairtools} stats {output.sort_pair} -o {output.pair_stat} ;\n"
		"echo {params.pairtools} stats {output.sort_pair} -o {output.pair_stat}"
		"    >> {params.logfile} ;\n\n"

rule figure_merge:
	input:
		statfile = absjoin("02.reads", "{dataset}", "{dataset}_clean", "quality_stat.txt"),
	output:
		pngfile = absjoin("02.reads", "{dataset}", "{dataset}_clean", 
						"{dataset}.sequencing_quality.png"),
		pdffile = absjoin("02.reads", "{dataset}", "{dataset}_clean", 
						"{dataset}.sequencing_quality.pdf"),
	params:
		python = software.get("python3"),
		figure_util = absjoin(scripts_dir, "figure_util.py"),
		base_content1 = absjoin("02.reads", "{dataset}", "{dataset}_clean_fastqc",
				"{dataset}_1P_fastqc", "Images", "per_base_sequence_content.png"),
		base_content2 = absjoin("02.reads", "{dataset}", "{dataset}_clean_fastqc",
				"{dataset}_2P_fastqc", "Images", "per_base_sequence_content.png"),
		GC_content1 = absjoin("02.reads", "{dataset}", "{dataset}_clean_fastqc",
				"{dataset}_1P_fastqc", "Images", "per_sequence_gc_content.png"),
		GC_content2 = absjoin("02.reads", "{dataset}", "{dataset}_clean_fastqc",
				"{dataset}_2P_fastqc", "Images", "per_sequence_gc_content.png"),
		base_quality1 = absjoin("02.reads", "{dataset}", "{dataset}_clean_fastqc",
				"{dataset}_1P_fastqc", "Images", "per_base_quality.png"),
		base_quality2 = absjoin("02.reads", "{dataset}", "{dataset}_clean_fastqc",
				"{dataset}_2P_fastqc", "Images", "per_base_quality.png"),
	shell:
		"{params.python} {params.figure_util} combine_figure --out-image {output.pngfile}"
		"    --ncol 3 {params.base_quality1} {params.base_content1} {params.GC_content1}"
		"    {params.base_quality2} {params.base_content2} {params.GC_content2} ;\n"
		
		"{params.python} {params.figure_util} combine_figure --out-image {output.pdffile}"
		"     --ncol 3 {params.base_quality1} {params.base_content1} {params.GC_content1}"
		"     {params.base_quality2} {params.base_content2} {params.GC_content2} ;\n"

rule read_stat:
	input:
		read1 = absjoin("02.reads", "{dataset}", "{dataset}_1P.fastq.gz"),
		read2 = absjoin("02.reads", "{dataset}", "{dataset}_2P.fastq.gz"),
	output:
		statfile = absjoin("02.reads", "{dataset}", "{dataset}_clean", "quality_stat.txt"),
	params:
		fastqc = software.get("fastqc"),
		mypython = software.get("python3"),
		r1_fastqc = "{dataset}_1P_fastqc",
		r2_fastqc = "{dataset}_2P_fastqc",
		seqkit = absjoin(src_dir, "seqkit"),
		clean_dir = absjoin("02.reads", "{dataset}", "{dataset}_clean"),
		raw_dir = absjoin("02.reads", "{dataset}", "{dataset}_clean_fastqc"),
	threads: 10
	shell:
		"{params.fastqc} --outdir {params.raw_dir} --extract"
		"    --threads {threads} {input.read1} {input.read2} ;\n"

		"cp {params.raw_dir}/{params.r1_fastqc}/Images/per_base_sequence_content.png"
		"    {params.clean_dir}/per_base_sequence_content_1.png ;\n"

		"cp {params.raw_dir}/{params.r2_fastqc}/Images/per_base_sequence_content.png"
		"    {params.clean_dir}/per_base_sequence_content_2.png ;\n"

		"cp {params.raw_dir}/{params.r1_fastqc}/Images/per_sequence_gc_content.png"
		"    {params.clean_dir}/per_sequence_GC_content_1.png ;\n"

		"cp {params.raw_dir}/{params.r2_fastqc}/Images/per_sequence_gc_content.png"
		"    {params.clean_dir}/per_sequence_GC_content_2.png ;\n"

		"cp {params.raw_dir}/{params.r1_fastqc}/Images/per_base_quality.png"
		"    {params.clean_dir}/per_base_sequence_quality_1.png ;\n"

		"cp {params.raw_dir}/{params.r2_fastqc}/Images/per_base_quality.png"
		"    {params.clean_dir}/per_base_sequence_quality_2.png ;\n"

		"cp {params.raw_dir}/{params.r1_fastqc}/Images/sequence_length_distribution.png"
		"    {params.clean_dir}/sequence_length_distribution_1.png ;\n"

		"cp {params.raw_dir}/{params.r2_fastqc}/Images/sequence_length_distribution.png"
		"    {params.clean_dir}/sequence_length_distribution_2.png ;\n"

		"""{params.seqkit} stat --all --basename --threads {threads} {input.read1} {input.read2} \
			  | sed 's/P.fastq.gz//g;s/P.fq.gz//g' | sed '1d' \
			  | awk -v OFS='\t' 'BEGIN{{print "SampleName","ReadNum","BaseCount(bp)", \
			  "ReadLength(bp)","Q20","Q30"}}{{print $1, $4, $5, $7, $14, $15}}' \
			  > {output.statfile} ;\n"""

rule read_clean:
	input:
		lambda wildcards: reads_dict.get(wildcards.dataset)
	output:
		read1 = absjoin("02.reads", "{dataset}", "{dataset}_1P.fastq.gz"),
		read2 = absjoin("02.reads", "{dataset}", "{dataset}_2P.fastq.gz"),
	params:
		java = software.get("java"),
		trimmomatic = software.get("trimmomatic"),
		adapter = get_adapter(config.get("platform")),
		baseout = absjoin("02.reads", "{dataset}", "{dataset}"),
		logfile = absjoin("02.reads", "work.sh"),
		trim_stat = absjoin("02.reads", "{dataset}", "{dataset}.trim.stats"),
	threads: 10
	shell:
		"{params.java} -Xmx35g -Djava.io.tmpdir=$PWD -jar {params.trimmomatic} PE"
		"    -threads {threads} -summary {params.trim_stat} -validatePairs {input}"
		"    ILLUMINACLIP:{params.adapter}:2:30:10:8:true LEADING:3 TRAILING:3"
		"    SLIDINGWINDOW:4:15 MINLEN:36 -baseout {params.baseout}.fastq.gz ;\n"
		"echo {params.java} -Xmx35g -Djava.io.tmpdir=$PWD -jar {params.trimmomatic} PE"
		"    -threads {threads} -summary {params.trim_stat} -validatePairs {input}"
		"    ILLUMINACLIP:{params.adapter}:2:30:10:8:true LEADING:3 TRAILING:3"
		"    SLIDINGWINDOW:4:15 MINLEN:36 -baseout {params.baseout}.fastq.gz >> {params.logfile} ;\n"

		"rm 02.reads/{wildcards.dataset}/*U.fastq.gz ;\n"

rule run_nr:
	input:
		protein = config.get("mrna_fa"),
		gene2mark = config.get("gene_mrna"),
	output:
		gene_nr = absjoin("01.ref", "NR", 'genome.gene.nraln'),
		nr_annot = absjoin("01.ref", "NR", "genome.go_annot"),
		wego = absjoin("01.ref", "NR", "genome.wego"),
		acc2go_output = absjoin("01.ref", "NR", "genome.ACC2GO"),
		classification = absjoin("01.ref", "NR", "genome.GO_classification_all.xls"),
	params:
		gotab = database.get("gotab"),
		perl = software.get("perl"),
		python = software.get("python3"),
		nr_pfix = absjoin("01.ref", "NR", "genome"),
		go_classify = software.get("go_classify"),
		nr2go = absjoin(scripts_dir, "nr2go.py"),
		access2go = database.get("accession2go"),
		modify = absjoin(dirname(software.get("go_classify")), "get_GO_anno.pl"), 
	threads: 24
	shell:
		"module load diamond/2.0.7 ;\n"
		"diamond blastx --db {nr_db} --query {input.protein} --threads {threads} -c1"
		"     --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart"
		"     send evalue bitscore stitle --more-sensitive --evalue 1e-5 --max-target-seqs 30"
		"     --salltitles --out {output.gene_nr} ;\n"

		"{params.python} {params.nr2go} --nr-align {output.gene_nr}"
		"     --access2go {params.access2go} --gotab {params.gotab}"
		"     --gene2mark {input.gene2mark} --outpfix {params.nr_pfix} ;\n"
	
		"{params.perl} {params.go_classify} {output.wego} {params.nr_pfix} ;\n"
		
		"{params.perl} {params.modify} -go {output.classification} -o {output.acc2go_output} ;\n"

rule run_kegg:
	input:
		protein = config.get("mrna_fa"),
		gene2mark = config.get("gene_mrna"),
	output:
		gene_ko = absjoin("01.ref", "KEGG", "genome.gene.ko"),
		gene_annot = absjoin("01.ref", "KEGG", "genome.gene.annot"),
		kegg2go_output = absjoin("01.ref", "KEGG", "genome.ACC2KEGG"),
	params:
		kobas_dir = software.get("kobas_dir"),
		kobas_annot = software.get("kobas_annot"),
		kobas_sqlite = software.get("kobas_sqlite"),
		python = software.get("python3"),
		kegg2go = absjoin(scripts_dir, "kegg2go.py"),
		kokeg = database.get("kokeg"),
		kegg_pfix = absjoin("01.ref", "KEGG", "genome"),
	threads: 24
	shell:
		"module load diamond/2.0.7 ;\n"
		"diamond blastx --db {kegg_db} --query {input.protein} -c1"
		"     --threads {threads} --outfmt 6 --evalue 1e-5 --more-sensitive"
		"     --max-target-seqs 15 --out {output.gene_ko} ;\n"

		"mypath=$PATH ;\n"
		"export PATH=/public/frasergen/RNA/pipeline/software/Anaconda2/anaconda2/bin:$PATH ;\n"
		"python {params.kobas_annot} -i {output.gene_ko} -t blastout:tab -s ko"
		"     -k {params.kobas_dir} -q {params.kobas_sqlite} -o {output.gene_annot} ;\n"
		"export PATH=/public/frasergen/PUB/software/python/Python-2.7.9/bin:$mypath ;\n"

		"{params.python} {params.kegg2go} --kokeg {params.kokeg}"
		"     --gene2mark {input.gene2mark} --kegg-annot {output.gene_annot}"
		"     --outpfix {params.kegg_pfix} ;\n"

#The output of this is directly used as the file to be provided.
#rule genome_fimo:
#	input:
#		inref = absjoin("01.ref", "genome.fa"),
#	output:
#		genome_motif = absjoin("01.ref", "motif_scan", "fimo.tsv"),
#	params:
#		fimo = software.get("fimo"),
#		python3 = software.get("python3"),
#		motif_db = config.get("jaspar_db"),
#		fimo_outdir = absjoin("01.ref", "motif_scan"),
#		parallel_fimo = absjoin(scripts_dir, "parallel_fimo.py"),
#	shell:
#		"{params.python3} {params.parallel_fimo} --infasta {input.inref} --meme_db {params.motif_db}"
#		"     --outdir {params.fimo_outdir} --fimo {params.fimo} ;\n"

rule ref_prepare:
	input:
		reference = absjoin(config.get("reference")),
		chrom_map = absjoin(config.get("chrom_map")),
	output:
		inref = absjoin("01.ref", "genome.fa"),
		infai = absjoin("01.ref", "genome.fa.fai"),
		reflen = absjoin("01.ref", "genome.chrsize"),
		genes = absjoin("01.ref", "genome_gene.bed"),
	params: 
		genome_gene = absjoin(config.get("gene_bed")),
		bwa = software.get("bwa"),
		python = software.get("python3"),
		fasize = software.get("fasize"),
		samtools = software.get("samtools"),
		chrom_trans = absjoin(scripts_dir, "chrom_trans.py"),
		logfile = absjoin("01.ref", "work.sh"),
	shell:
		"{params.python} {params.chrom_trans} --inref {input.reference}"
		"    --gene-bed {params.genome_gene} --chrom-map {input.chrom_map}"
		"    --outfa {output.inref} --outgene {output.genes} ;\n"
		"echo {params.python} {params.chrom_trans} --inref {input.reference}"
		"    --gene-bed {params.genome_gene} --chrom-map {input.chrom_map}"
		"    --outfa {output.inref} --outgene {output.genes} > {params.logfile} ;\n\n"
		#"ln -s {input.reference} {output.inref} ;\n"
		#"ln -s {input.genome_gene} {output.genes} ;\n"

		"{params.samtools} faidx {output.inref} ;\n"
		"echo {params.samtools} faidx {output.inref} >> {params.logfile} ;\n\n"

		"{params.bwa} index {output.inref} ;\n"
		"echo {params.bwa} index {output.inref} >> {params.logfile} ;\n\n"
			
		"{params.fasize} -detailed -tab {output.inref} > {output.reflen} ;\n"
		"echo {params.fasize} -detailed -tab {output.inref} \> {output.reflen} >> {params.logfile} ;\n\n"

