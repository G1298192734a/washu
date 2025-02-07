import sys, snakemake, yaml
import os, time, pysam, math
from pathlib import Path
absjoin = lambda *intuple: os.path.abspath(os.path.join(*intuple))

snake_dir = Path("/public/frasergen/3D/pipeline/Interactome/hic3_workflow/hic3_biorep.py").resolve().parent
inyaml = snake_dir.joinpath("software.yaml")
biorep_dir = snake_dir.joinpath("biorep")
src_dir = snake_dir.joinpath("src")
scripts_dir = snake_dir.joinpath("scripts")
software = yaml.safe_load(open(inyaml))

spl_minus = config.get("sample_minus").replace(" ", "").split(",")
single_dir = Path(config.get("single_sample_dir"))
genome_res = config.get("genome_res")
compartment_res = config.get("compartment_res")
tad_res = config.get("tad_res")
loop_res = config.get("loop_res")
sample_order = config.get("sample_order").replace(" ", "").split(",")

dir1 = Path("01.diff_heatmap").resolve()
dir2 = Path("02.diff_compartment").resolve()
dir3 = Path("03.diff_TAD").resolve()
dir4 = Path("04.diff_loop").resolve()
snakemake.utils.makedirs([dir1.joinpath("01.genome_diffhm"), 
						dir1.joinpath("02.bychr_diffhm"),  
						dir1.joinpath("03.diff_pscurve"), 
						dir2.joinpath("All"), dir3, dir4, 
						dir4.joinpath("01.mustache", "All"), 
						dir4.joinpath("02.fithic", "All")])

def get_switch_name(wildcards):
	case, control = wildcards.spl_pair.split(".minus.")
	return dir2.joinpath(wildcards.spl_pair, f"{control}_to_{case}.switch.xls")

rule all:
	input:
		expand(dir1.joinpath("01.genome_diffhm", "{spl_pair}", "{spl_pair}.genome.png"), 
				spl_pair=spl_minus),
#		dir1.joinpath("03.diff_pscurve", "decay_curve.png"),
#		dir2.joinpath("All", "saddle_plot_all.png"),
#		expand(dir2.joinpath("{spl_pair}", "{spl_pair}.switch.xls"), spl_pair=spl_minus),
#		expand(dir2.joinpath("{spl_pair}", "{spl_pair}.A2B", "{spl_pair}.A2B.enrichGO.png"),
#						spl_pair=spl_minus),
#		expand(dir2.joinpath("{spl_pair}", "switch_plot.{compartment_res}", "{spl_pair}.chr1.pdf"), 
#						spl_pair=spl_minus, compartment_res=compartment_res),
#		expand(dir3.joinpath("02.hicfindtad", "{spl_pair}", "boundary_enrich.png"),
#						spl_pair=spl_minus),
#		expand(dir3.joinpath("02.hicfindtad", "{spl_pair}", "{spl_pair}.classify.xls"), 
#						spl_pair=spl_minus),
#		expand(dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.domain_ata.png"),
#						spl_pair=spl_minus),
#		expand(dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.boundary_ata.png"),
#						spl_pair=spl_minus),
#		expand(dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.classify.xls"),
#						spl_pair=spl_minus),
##		expand(dir4.joinpath("02.fithic", "{spl_pair}", "{spl_pair}-cis.venn.png"), 
##						spl_pair=spl_minus),
#		expand(dir4.joinpath("01.mustache", "{spl_pair}", "{spl_pair}-share_loops.xls"),
#						spl_pair=spl_minus),
#		expand(dir4.joinpath("01.mustache", "{spl_pair}", "{spl_pair}.APA.png"),
#						spl_pair=spl_minus),
#		dir4.joinpath("01.mustache", "All", "loop_length_all.png"),
#		dir4.joinpath("02.fithic", "All", "loop_length_all.png"),

rule all_loop_fithic:
	input: 
		loops = expand(single_dir.joinpath("07.loops", "03.fithic", "{spl}.{loop_res}", "04.circos", 
						"{spl}_{loop_res}.cis_interact.bedpe"), spl=sample_order, loop_res=loop_res),
	output:
		loop_length_png = dir4.joinpath("02.fithic", "All", "loop_length_all.png"),
		loop_length_pdf = dir4.joinpath("02.fithic", "All", "loop_length_all.pdf"),
	params:
		rscript = software.get("Rscript"),
		outpfix = dir4.joinpath("02.fithic", "All", "loop_length_all"),
		loop_len_script = biorep_dir.joinpath("loop_length_all.r"),
		sample_names = config.get("sample_order").replace(" ", ""),
		loop_bedpe = ",".join([str(single_dir.joinpath("07.loops", "03.fithic", f"{spl}.{loop_res}", 
								"04.circos", f"{spl}_{loop_res}.cis_interact.bedpe")) for spl in sample_order]),
	shell:
		"{params.rscript} {params.loop_len_script} --names {params.sample_names} "
		"    --loop_bedpe {params.loop_bedpe} --outpfix {params.outpfix} ;\n"

rule diff_fithic:
	input:  
		acc2go = single_dir.joinpath("01.ref", "NR", "genome.ACC2GO"),
		kegg2go = single_dir.joinpath("01.ref", "KEGG", "genome.ACC2KEGG"),
		gene_bed = single_dir.joinpath("01.ref", "genome_gene.bed"),
		cisfiles = lambda wildcards: single_dir.joinpath("07.loops", "03.fithic", 
							wildcards.spl_pair.split(".minus.")[0]+f".{loop_res}", "04.circos",
							wildcards.spl_pair.split(".minus.")[0]+f"_{loop_res}.cis.interact.gz"),
		loop_cool = lambda wildcards: [single_dir.joinpath("07.loops", "01.prepare", 
							f"{spl}_{loop_res}.cool") for spl in wildcards.spl_pair.split(".minus.")],
	output:
		cis_venn = dir4.joinpath("02.fithic", "{spl_pair}", "{spl_pair}-cis.venn.png"),
	params:
		python3 = software.get("python3"),
		rscript = software.get("Rscript"),
		bedtools = software.get("bedtools"),
		enrich = biorep_dir.joinpath("enrichment.R"),
		fithic_compare = biorep_dir.joinpath("fithic_compare.py"),
		cluster_profiler = biorep_dir.joinpath("clusterProfiler.R"),
		splnames = lambda wildcards: wildcards.spl_pair.replace(".minus.", ","),
		control_name = lambda wildcards: wildcards.spl_pair.split(".minus.")[1],
		treat_name = lambda wildcards: wildcards.spl_pair.split(".minus.")[0],
		outpfix = lambda wildcards: dir4.joinpath("02.fithic", wildcards.spl_pair, wildcards.spl_pair),
		cis_interact = lambda wildcards: ",".join([str(single_dir.joinpath("07.loops", "03.fithic", 
								f"{spl}.{loop_res}", "04.circos", f"{spl}_{loop_res}.cis.interact.gz")) 
								for spl in wildcards.spl_pair.split(".minus.")]),
		trans_interact = lambda wildcards: ",".join([str(single_dir.joinpath("07.loops", "03.fithic", 
								f"{spl}.{loop_res}", "04.circos", f"{spl}_{loop_res}.trans.interact.gz"))
								for spl in wildcards.spl_pair.split(".minus.")]),
		treat_specific = lambda wildcards: dir4.joinpath("02.fithic", wildcards.spl_pair, 
						wildcards.spl_pair+"-"+wildcards.spl_pair.split(".minus.")[0] + "_specific"),
		control_specific = lambda wildcards: dir4.joinpath("02.fithic", wildcards.spl_pair,
						wildcards.spl_pair+"-"+wildcards.spl_pair.split(".minus.")[1] + "_specific"),
		share = lambda wildcards: dir4.joinpath("02.fithic", wildcards.spl_pair,
						wildcards.spl_pair + "-share"),
	shell:
		"{params.python3} {params.fithic_compare} --cis {params.cis_interact}"
		"     --trans {params.trans_interact} --sample {params.splnames}"
		"     --res {loop_res} --output {params.outpfix} ;\n"

		"{params.bedtools} pairtobed -a {params.treat_specific}_cis.bedpe -b {input.gene_bed} |"
		"     cut -f 10 | sort -u > {params.treat_specific}_cis.genelist ;\n"
	
		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"     --genelist_file {params.treat_specific}_cis.genelist -o {params.treat_specific}_cis"
		"     -p {wildcards.spl_pair}-{params.treat_name}_cis ;\n"
		
		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {params.treat_specific}_cis.genelist"
		"     {params.treat_specific}_cis/{wildcards.spl_pair}-{params.treat_name}_cis ;\n"

		"{params.bedtools} pairtobed -a {params.control_specific}_cis.bedpe -b {input.gene_bed} |"
		"     cut -f 10 | sort -u > {params.control_specific}_cis.genelist ;\n"
		
		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"     --genelist_file {params.control_specific}_cis.genelist -o {params.control_specific}_cis"
		"     -p {wildcards.spl_pair}-{params.control_name}_cis ;\n"
	
		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {params.control_specific}_cis.genelist"
		"     {params.control_specific}_cis/{wildcards.spl_pair}-{params.control_name}_cis ;\n"

		"{params.bedtools} pairtobed -a {params.share}_cis.bedpe -b {input.gene_bed} |"
		"     cut -f 10 | sort -u > {params.share}_cis.genelist ;\n"

		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"     --genelist_file {params.share}_cis.genelist -o {params.share}_cis"
		"     -p {wildcards.spl_pair}-share_cis ;\n"

		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {params.share}_cis.genelist"
		"     {params.share}_cis/{wildcards.spl_pair}-share_cis ;\n"

		"{params.bedtools} pairtobed -a {params.treat_specific}_trans.bedpe -b {input.gene_bed} |"
		"     cut -f 10 | sort -u > {params.treat_specific}_trans.genelist ;\n"

		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"     --genelist_file {params.treat_specific}_trans.genelist -o {params.treat_specific}_trans"
		"     -p {wildcards.spl_pair}-{params.treat_name}_trans ;\n"

		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {params.treat_specific}_trans.genelist"
		"     {params.treat_specific}_trans/{wildcards.spl_pair}-{params.treat_name}_trans ;\n"

		"{params.bedtools} pairtobed -a {params.control_specific}_trans.bedpe -b {input.gene_bed} |"
		"     cut -f 10 | sort -u > {params.control_specific}_trans.genelist ;\n"

		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"     --genelist_file {params.control_specific}_trans.genelist -o {params.control_specific}_trans"
		"     -p {wildcards.spl_pair}-{params.control_name}_trans ;\n"
	
		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {params.control_specific}_trans.genelist"
		"     {params.control_specific}_trans/{wildcards.spl_pair}-{params.control_name}_trans ;\n"

		"{params.bedtools} pairtobed -a {params.share}_trans.bedpe -b {input.gene_bed} |"
		"     cut -f 10 | sort -u > {params.share}_trans.genelist ;\n"
		
		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"     --genelist_file {params.share}_trans.genelist -o {params.share}_trans"
		"     -p {wildcards.spl_pair}-share_trans ;\n"

		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {params.share}_trans.genelist"
		"     {params.share}_trans/{wildcards.spl_pair}-share_trans ;\n"

def loop_path(spl_pair, control_treat, single_dir):
	treat, control = spl_pair.split(".minus.")
	if control_treat == "control":
		return single_dir.joinpath("07.loops", "02.mustache", control, f"{control}.combine.mustache")
	else:
		return single_dir.joinpath("07.loops", "02.mustache", treat, f"{treat}.combine.mustache")

rule all_loop_mustache:
	input:
		loops = expand(single_dir.joinpath("07.loops", "02.mustache", "{spl}", "{spl}.combine.bedpe"), 
						spl=sample_order),
	output:
		loop_length_png = dir4.joinpath("01.mustache", "All", "loop_length_all.png"),
		loop_length_pdf = dir4.joinpath("01.mustache", "All", "loop_length_all.pdf"),
	params:
		rscript = software.get("Rscript"),
		outpfix = dir4.joinpath("01.mustache", "All", "loop_length_all"),
		loop_len_script = biorep_dir.joinpath("loop_length_all.r"),
		sample_names = config.get("sample_order").replace(" ", ""),
		loop_bedpe = ",".join([str(single_dir.joinpath("07.loops", "02.mustache", spl, 
								f"{spl}.combine.bedpe")) for spl in sample_order]),
	shell:
		"{params.rscript} {params.loop_len_script} --names {params.sample_names} "
		"    --loop_bedpe {params.loop_bedpe} --outpfix {params.outpfix} ;\n"

rule diff_mustache:
	input:
		gene_bed = single_dir.joinpath("01.ref", "genome_gene.bed"),
		acc2go = single_dir.joinpath("01.ref", "NR", "genome.ACC2GO"),
		kegg2go = single_dir.joinpath("01.ref", "KEGG", "genome.ACC2KEGG"),
		control_loops = lambda wildcards: loop_path(wildcards.spl_pair, "control", single_dir),
		treat_loops = lambda wildcards: loop_path(wildcards.spl_pair, "treat", single_dir),
		loop_cool = lambda wildcards: [single_dir.joinpath("07.loops", "01.prepare", 
						f"{spl}_{loop_res}.cool") for spl in wildcards.spl_pair.split(".minus.")],
	output:
		mustache_share = dir4.joinpath("01.mustache", "{spl_pair}", "{spl_pair}-share_loops.xls"),
		control_specific = temp(dir4.joinpath("01.mustache", "{spl_pair}", "control_specific_loops.bedpe")),
		treat_specific = temp(dir4.joinpath("01.mustache", "{spl_pair}", "treat_specific_loops.bedpe")),
		share_loops = temp(dir4.joinpath("01.mustache", "{spl_pair}", "share_loops.bedpe")),
		cloop_chm = temp(dir4.joinpath("01.mustache", "{spl_pair}", "control_loop_control_hm.apa")), 
		cloop_thm = temp(dir4.joinpath("01.mustache", "{spl_pair}", "control_loop_treat_hm.apa")),
		tloop_chm = temp(dir4.joinpath("01.mustache", "{spl_pair}", "treat_loop_control_hm.apa")),
		tloop_thm = temp(dir4.joinpath("01.mustache", "{spl_pair}", "treat_loop_treat_hm.apa")),
		sloop_chm = temp(dir4.joinpath("01.mustache", "{spl_pair}", "share_loop_control_hm.apa")), 
		sloop_thm = temp(dir4.joinpath("01.mustache", "{spl_pair}", "share_loop_treat_hm.apa")),
		apa_png = dir4.joinpath("01.mustache", "{spl_pair}", "{spl_pair}.APA.png"),
		apa_pdf = dir4.joinpath("01.mustache", "{spl_pair}", "{spl_pair}.APA.pdf"),
	params:
		python3 = software.get("python3"),
		rscript = software.get("Rscript"),
		coolpup = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/coolpup.py",
		bedtools = software.get("bedtools"),
		enrich = biorep_dir.joinpath("enrichment.R"),
		plotpup = src_dir.joinpath("plotpup", "plotpup.py"),
		mustache_compare = biorep_dir.joinpath("mustache_compare.py"),
		cluster_profiler = biorep_dir.joinpath("clusterProfiler.R"),
		control_name = lambda wildcards: wildcards.spl_pair.split(".minus.")[1],
		treat_name = lambda wildcards: wildcards.spl_pair.split(".minus.")[0],
		control_specific = lambda wildcards: dir4.joinpath("01.mustache", wildcards.spl_pair, 
					wildcards.spl_pair+"-"+wildcards.spl_pair.split(".minus.")[1] + "_specific"),
		treat_specific = lambda wildcards: dir4.joinpath("01.mustache", wildcards.spl_pair, 
					wildcards.spl_pair+"-"+wildcards.spl_pair.split(".minus.")[0] + "_specific"),	
		outdir = lambda wildcards: dir4.joinpath("01.mustache", wildcards.spl_pair),
	shell:
		"{params.python3} {params.mustache_compare} --control_loops {input.control_loops}"
		"     --treat_loops {input.treat_loops} --control_name {params.control_name}"
		"     --treat_name {params.treat_name} --outdir {params.outdir} ;\n"

		"{params.coolpup} {input.loop_cool[1]} {output.control_specific} --outname {output.cloop_chm} ;\n"	

		"{params.coolpup} {input.loop_cool[1]} {output.treat_specific} --outname {output.tloop_chm} ;\n"

		"{params.coolpup} {input.loop_cool[1]} {output.share_loops} --outname {output.sloop_chm} ;\n"

		"{params.coolpup} {input.loop_cool[0]} {output.control_specific} --outname {output.cloop_thm} ;\n"

		"{params.coolpup} {input.loop_cool[0]} {output.treat_specific} --outname {output.tloop_thm} ;\n"

		"{params.coolpup} {input.loop_cool[0]} {output.share_loops} --outname {output.sloop_thm} ;\n"

		"{params.plotpup} {output.cloop_chm} {output.cloop_thm} {output.tloop_chm} {output.tloop_thm}"
		"     {output.sloop_chm} {output.sloop_thm} --output {output.apa_png}"
		"     --dpi 600 --n_cols 2 --col_names {params.control_name},{params.treat_name}"
		"     --row_names {params.control_name}_specific,{params.treat_name}_specific,Shared ;\n"

		"{params.plotpup} {output.cloop_chm} {output.cloop_thm} {output.tloop_chm} {output.tloop_thm}"
		"     {output.sloop_chm} {output.sloop_thm} --output {output.apa_pdf}"
		"     --dpi 600 --n_cols 2 --col_names {params.control_name},{params.treat_name}"
		"     --row_names {params.control_name}_specific,{params.treat_name}_specific,Shared ;\n"

		"sed '1d' {params.control_specific}_loops.xls | cut -f 1-6 |"
		"     {params.bedtools} pairtobed -a stdin -b {input.gene_bed} |"
		"     cut -f 10 | sort -u > {params.control_specific}.genelist ;\n"

		"sed '1d' {params.treat_specific}_loops.xls | cut -f 1-6 |"
		"     {params.bedtools} pairtobed -a stdin -b {input.gene_bed} |"
		"     cut -f 10 | sort -u > {params.treat_specific}.genelist ;\n"
		
		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"     --genelist_file {params.control_specific}.genelist -o {params.control_specific}"
		"     -p {wildcards.spl_pair}-{params.control_name} ;\n"
		
		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {params.control_specific}.genelist"
		"     {params.control_specific}/{wildcards.spl_pair}-{params.control_name} ;\n"

		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"     --genelist_file {params.treat_specific}.genelist -o {params.treat_specific}"
		"     -p {wildcards.spl_pair}-{params.treat_name} ;\n"
	
		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {params.treat_specific}.genelist"
		"     {params.treat_specific}/{wildcards.spl_pair}-{params.treat_name} ;\n"

rule diff_diamond_domain:
	input:
		tad_basic_cool = lambda wildcards: [single_dir.joinpath("06.tad", "01.prepare",
						f"{spl}.basic.cool") for spl in wildcards.spl_pair.split(".minus.")],
		tad_region = lambda wildcards: expand(single_dir.joinpath("06.tad", "02.cooltools", 
						"{spl}.{tad_res}", "TAD_genome.bed"), 
						spl=wildcards.spl_pair.split(".minus."), tad_res=tad_res),
	output:
		diamond_classify = dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.classify.xls"),
		classify_stat_png = dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.classify.stat.png"),
		stable_region = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Stable.bed")),
		merge_region = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Merge.bed")),
		split_region = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Split.bed")),
		rearr_region = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Rearrangement.bed")),
		stable_control = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Stable_control.ata")),
		merge_control = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Merge_control.ata")),
		split_control = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Split_control.ata")),
		rearr_control = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Rearrangement_control.ata")),
		stable_treat = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Stable_treat.ata")),
		merge_treat = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Merge_treat.ata")),
		split_treat = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Split_treat.ata")),
		rearr_treat = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "Rearrangement_treat.ata")),
		domain_ata_png = dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.domain_ata.png"),
		domain_ata_pdf = dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.domain_ata.pdf"),
	params:
		python3 = software.get("python3"),
		coolpup = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/coolpup.py",
		tad_classify = biorep_dir.joinpath("tad_classify.py"),
		tad_classify_size = biorep_dir.joinpath("tad_classify_size.py"),
		name_pair = lambda wildcards: wildcards.spl_pair.split(".minus."),
		classify_outpfix = lambda wildcards: dir3.joinpath("01.cooltools", 
						wildcards.spl_pair, wildcards.spl_pair+".classify"),
		plotpup = src_dir.joinpath("plotpup", "plotpup.py"),
	shell:
		"{params.python3} {params.tad_classify} --treat-tad {input.tad_region[0]}"
		"    --control-tad {input.tad_region[1]} --treat-name {params.name_pair[0]}"
		"    --control-name {params.name_pair[1]} --outpfix {params.classify_outpfix} ;\n"

		"{params.python3} {params.tad_classify_size} --tad-classify {output.diamond_classify}"
		"    --outpfix {params.classify_outpfix} ;\n"

		"{params.coolpup} --local --rescale --minsize 60000 --maxsize 3000000 --n_proc 32"
		"    --outname {output.stable_control} --seed 100 {input.tad_basic_cool[1]} {output.stable_region} ;\n"

		"{params.coolpup} --local --rescale --minsize 60000 --maxsize 3000000 --n_proc 32"
		"    --outname {output.merge_control} --seed 100 {input.tad_basic_cool[1]} {output.merge_region} ;\n"

		"{params.coolpup} --local --rescale --minsize 60000 --maxsize 3000000 --n_proc 32"
		"    --outname {output.split_control} --seed 100 {input.tad_basic_cool[1]} {output.split_region} ;\n"

		"{params.coolpup} --local --rescale --minsize 60000 --maxsize 3000000 --n_proc 32"
		"    --outname {output.rearr_control} --seed 100 {input.tad_basic_cool[1]} {output.rearr_region} ;\n"
		
		"{params.coolpup} --local --rescale --minsize 60000 --maxsize 3000000 --n_proc 32"
		"    --outname {output.stable_treat} --seed 100 {input.tad_basic_cool[0]} {output.stable_region} ;\n"
		
		"{params.coolpup} --local --rescale --minsize 60000 --maxsize 3000000 --n_proc 32"
		"    --outname {output.merge_treat} --seed 100 {input.tad_basic_cool[0]} {output.merge_region} ;\n"
		
		"{params.coolpup} --local --rescale --minsize 60000 --maxsize 3000000 --n_proc 32"
		"    --outname {output.split_treat} --seed 100 {input.tad_basic_cool[0]} {output.split_region} ;\n"

		"{params.coolpup} --local --rescale --minsize 60000 --maxsize 3000000 --n_proc 32"
		"    --outname {output.rearr_treat} --seed 100 {input.tad_basic_cool[0]} {output.rearr_region} ;\n"

		"{params.plotpup} --scale log --output {output.domain_ata_png} --dpi 600"
		"    --row_names Stable,Merge,Split,Rearrangement --col_names {params.name_pair[1]},{params.name_pair[0]}"
		"    --n_cols 2 {output.stable_control} {output.stable_treat} {output.merge_control} {output.merge_treat}"
		"    {output.split_control} {output.split_treat} {output.rearr_control} {output.rearr_treat} ;\n"

		"{params.plotpup} --scale log --output {output.domain_ata_pdf}"
		"    --row_names Stable,Merge,Split,Rearrangement --col_names {params.name_pair[1]},{params.name_pair[0]}"
		"    --n_cols 2 {output.stable_control} {output.stable_treat} {output.merge_control} {output.merge_treat}"
		"    {output.split_control} {output.split_treat} {output.rearr_control} {output.rearr_treat} ;\n"

rule diff_diamond_boundary:
	input:
		gene_bed = single_dir.joinpath("01.ref", "genome_gene.bed"),
		acc2go = single_dir.joinpath("01.ref", "NR", "genome.ACC2GO"),
		kegg2go = single_dir.joinpath("01.ref", "KEGG", "genome.ACC2KEGG"),
		tad_basic_cool = lambda wildcards: [single_dir.joinpath("06.tad", "01.prepare", 
					f"{spl}.basic.cool") for spl in wildcards.spl_pair.split(".minus.")],
		diamond_boundary = lambda wildcards: [single_dir.joinpath("06.tad", "02.cooltools", 
					f"{spl}.{tad_res}", f"{spl}.{tad_res}_boundaries.bed")
					for spl in wildcards.spl_pair.split(".minus.")],
	output:
		venn_diagram = dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.boundary_venn.png"),
		share_boundaries = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "shared_boundaries.bed")),
		control_boundaries = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "control_specific_boundaries.bed")),
		treat_boundaries = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "treat_specific_boundaries.bed")),
		cc_ata = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "control_boundary.control_hm.ata")),
		ct_ata = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "control_boundary.treat_hm.ata")),
		tc_ata = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "treat_boundary.control_hm.ata")),
		tt_ata = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "treat_boundary.treat_hm.ata")),
		sc_ata = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "shared_boundary.treat_hm.ata")),
		st_ata = temp(dir3.joinpath("01.cooltools", "{spl_pair}", "shared_boundary.control_hm.ata")),
		all_boundary = dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.boundary_overlap.xls"),
		boundary_ata_png = dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.boundary_ata.png"),
		boundary_ata_pdf = dir3.joinpath("01.cooltools", "{spl_pair}", "{spl_pair}.boundary_ata.pdf"),
		control_genelist = dir3.joinpath("01.cooltools", "{spl_pair}", "control_specific_boundary.genelist"),
		treat_genelist = dir3.joinpath("01.cooltools", "{spl_pair}", "treat_specific_boundary.genelist"),
	params:
		python3 = software.get("python3"),
		coolpup = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/coolpup.py",
		bedtools = software.get("bedtools"),
		rscript = software.get("Rscript"),
		enrich = biorep_dir.joinpath("enrichment.R"),
		cluster_profiler = biorep_dir.joinpath("clusterProfiler.R"),
		plotpup = src_dir.joinpath("plotpup", "plotpup.py"),
		bATA_outdir = lambda wildcards: dir3.joinpath("01.cooltools", wildcards.spl_pair),
		diff_boundary = biorep_dir.joinpath("diamond_diff_boundary.py"),
		name_pair = lambda wildcards: wildcards.spl_pair.split(".minus."),
		enrich_outdir = lambda wildcards: [dir3.joinpath("01.cooltools", wildcards.spl_pair, f"{spl}_specific") 
						for spl in wildcards.spl_pair.split(".minus.")],
	threads: 32
	shell:
		"{params.python3} {params.diff_boundary} --control_name {params.name_pair[1]}"
		"    --treat_name {params.name_pair[0]} --control_boundary {input.diamond_boundary[1]}"
		"    --treat_boundary {input.diamond_boundary[0]} --outdir {params.bATA_outdir}"
		"    --gene_bed {input.gene_bed} ;\n"	

		"{params.coolpup} --local --n_proc {threads} --outname {output.cc_ata} --ignore_diags 0"
		"    --seed 100 --pad 500 {input.tad_basic_cool[1]} {output.control_boundaries} ;\n"

		"{params.coolpup} --local --n_proc {threads} --outname {output.ct_ata} --ignore_diags 0"
		"    --seed 100 --pad 500 {input.tad_basic_cool[0]} {output.control_boundaries} ;\n"

		"{params.coolpup} --local --n_proc {threads} --outname {output.tc_ata} --ignore_diags 0"
		"    --seed 100 --pad 500 {input.tad_basic_cool[1]} {output.treat_boundaries} ;\n"

		"{params.coolpup} --local --n_proc {threads} --outname {output.tt_ata} --ignore_diags 0"
		"    --seed 100 --pad 500 {input.tad_basic_cool[0]} {output.treat_boundaries} ;\n"

		"{params.coolpup} --local --n_proc {threads} --outname {output.sc_ata} --ignore_diags 0"
		"    --seed 100 --pad 500 {input.tad_basic_cool[1]} {output.share_boundaries} ;\n"

		"{params.coolpup} --local --n_proc {threads} --outname {output.st_ata} --ignore_diags 0"
		"    --seed 100 --pad 500 {input.tad_basic_cool[0]} {output.share_boundaries} ;\n"

		"{params.plotpup} --scale log --output {output.boundary_ata_png} --dpi 600 --n_cols 2"
		"    --col_names {params.name_pair[1]},{params.name_pair[0]}"
		"    --row_names {params.name_pair[1]}_specific,{params.name_pair[0]}_specific,shared"
		"    {output.cc_ata} {output.ct_ata} {output.tc_ata} {output.tt_ata} {output.sc_ata} {output.st_ata} ;\n"

		"{params.plotpup} --scale log --output {output.boundary_ata_pdf} --n_cols 2"
		"    --col_names {params.name_pair[1]},{params.name_pair[0]}"
		"    --row_names {params.name_pair[1]}_specific,{params.name_pair[0]}_specific,shared"
		"    {output.cc_ata} {output.ct_ata} {output.tc_ata} {output.tt_ata} {output.sc_ata} {output.st_ata} ;\n"
	
		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"    --genelist_file {output.control_genelist} -o {params.enrich_outdir[1]}"
		"    -p {params.name_pair[1]}_specific ;\n"

		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {output.control_genelist}"
		"    {params.enrich_outdir[1]}/{params.name_pair[1]}_specific ;\n"

		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"    --genelist_file {output.treat_genelist} -o {params.enrich_outdir[0]}"
		"    -p {params.name_pair[0]}_specific ;\n"

		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {output.treat_genelist}"
		"    {params.enrich_outdir[0]}/{params.name_pair[0]}_specific ;\n"
	
rule diff_hicfindtad:
	input:
		infai = single_dir.joinpath("01.ref", "genome.fa.fai"),
		chrsize = single_dir.joinpath("01.ref", "genome.chrsize"),
		in_boundary = lambda wildcards: [single_dir.joinpath("06.tad", "03.hicFindTAD", 
					f"{spl}.{tad_res}", f"{spl}.{tad_res}_boundaries.bed")
					for spl in wildcards.spl_pair.split(".minus.")],
		in_bw = lambda wildcards: [single_dir.joinpath("06.tad", "03.hicFindTAD",
					f"{spl}.{tad_res}", f"{spl}.{tad_res}.score.bw")
					for spl in wildcards.spl_pair.split(".minus.")],
		tad_cool = lambda wildcards: expand(single_dir.joinpath("06.tad", "01.prepare", 
					"{spl}_{tad_res}.tad.cool"), 
					spl=wildcards.spl_pair.split(".minus."), tad_res=tad_res),
		in_domain = lambda wildcards: expand(single_dir.joinpath("06.tad", "03.hicFindTAD", 
					"{spl}.{tad_res}", "{spl}.{tad_res}_domains.bed"), 
					spl=wildcards.spl_pair.split(".minus."), tad_res=tad_res),
		tad_region = lambda wildcards: expand(single_dir.joinpath("06.tad", "03.hicFindTAD", 
					"{spl}.{tad_res}", "TAD_genome.bed"), 
					spl=wildcards.spl_pair.split(".minus."), tad_res=tad_res),
	output:
		tad_figure = dir3.joinpath("02.hicfindtad", "{spl_pair}", "{spl_pair}.boundary.png"),
		classify_output = dir3.joinpath("02.hicfindtad", "{spl_pair}", "{spl_pair}.classify.xls"),
		enrich_figure = dir3.joinpath("02.hicfindtad", "{spl_pair}", "boundary_enrich.png"),
	params:
		python3 = software.get("python3"),
		rscript = software.get("Rscript"),
		tad_len = biorep_dir.joinpath("tad_length.r"),
		tad_classify = biorep_dir.joinpath("tad_classify.py"),
		tad_classify_size = biorep_dir.joinpath("tad_classify_size.py"),
		tad_enrich = biorep_dir.joinpath("tad_enrich.py"),
		name_pair = lambda wildcards: wildcards.spl_pair.split(".minus."),
		name_pair_str = lambda wildcards: wildcards.spl_pair.replace(".minus.", ","),
		diff_boundary = biorep_dir.joinpath("diff_boundary.py"),
		outdir = lambda wildcards: dir3.joinpath("02.hicfindtad", wildcards.spl_pair),
		outpfix = lambda wildcards: dir3.joinpath("02.hicfindtad", wildcards.spl_pair, wildcards.spl_pair),
		tad_bed = lambda wildcards: ",".join([str(single_dir.joinpath("06.tad", "03.hicFindTAD",
					f"{spl}.{tad_res}", "TAD_genome.bed")) 
					for spl in wildcards.spl_pair.split(".minus.")]),
		hicDifferentialTAD = Path(software.get("hicexplorer")).parent.joinpath("hicDifferentialTAD"),
	threads: 24
	shell:
		"{params.python3} {params.diff_boundary} --infai {input.infai}"
		"    --control-boundary {input.in_boundary[1]} --case-boundary {input.in_boundary[0]}"
		"    --control-name {params.name_pair[1]} --case-name {params.name_pair[0]}"
		"    --outpfix {params.outpfix} ;\n"

		"{params.rscript} {params.tad_len} --tad_beds {params.tad_bed}"
		"    --splnames {params.name_pair_str} --outdir {params.outdir} ;\n"

		"{params.hicDifferentialTAD} --targetMatrix {input.tad_cool[0]} --controlMatrix {input.tad_cool[1]}"
		"    --tadDomains {input.in_domain[0]} --outFileNamePrefix {params.outdir}/{params.name_pair[0]}"
		"    --threads {threads} ;\n"
		
		"{params.hicDifferentialTAD} --targetMatrix {input.tad_cool[0]} --controlMatrix {input.tad_cool[1]}"
		"    --tadDomains {input.in_domain[1]} --outFileNamePrefix {params.outdir}/{params.name_pair[1]}"
		"    --threads {threads} ;\n"

		"{params.python3} {params.tad_classify} --treat-tad {input.tad_region[0]}"
		"    --control-tad {input.tad_region[1]} --treat-name {params.name_pair[0]}"
		"    --control-name {params.name_pair[1]} --outpfix {params.outpfix}.classify ;\n"

		"{params.python3} {params.tad_classify_size} --tad-classify {output.classify_output}"
		"    --outpfix {params.outpfix}.classify ;\n"

		"{params.python3} {params.tad_enrich} --chrsize {input.chrsize}"
		"    --control_boundary {input.in_boundary[1]} --control_bw {input.in_bw[1]}"
		"    --treat_boundary {input.in_boundary[0]} --treat_bw {input.in_bw[0]}"
		"    --control_name {params.name_pair[1]} --treat_name {params.name_pair[0]}"
		"    --outdir {params.outdir} ;\n"

rule enrich_switch:
	input:
		gene_bed = single_dir.joinpath("01.ref", "genome_gene.bed"),
		acc2go = single_dir.joinpath("01.ref", "NR", "genome.ACC2GO"),
		kegg2go = single_dir.joinpath("01.ref", "KEGG", "genome.ACC2KEGG"),
		switch_table = dir2.joinpath("{spl_pair}", "{spl_pair}.switch.xls"),
	output:
		A2B_bed = temp(dir2.joinpath("{spl_pair}", "{spl_pair}.A2B.bed")),
		B2A_bed = temp(dir2.joinpath("{spl_pair}", "{spl_pair}.B2A.bed")),
		A2B_gene = dir2.joinpath("{spl_pair}", "{spl_pair}.A2B.genelist"),	
		B2A_gene = dir2.joinpath("{spl_pair}", "{spl_pair}.B2A.genelist"),	
		A2B_enrich_go = dir2.joinpath("{spl_pair}", "{spl_pair}.A2B", "{spl_pair}.A2B.enrichGO.png"),
		B2A_enrich_go = dir2.joinpath("{spl_pair}", "{spl_pair}.B2A", "{spl_pair}.B2A.enrichGO.png"),
	params:
		rscript = software.get("Rscript"),
		enrich = biorep_dir.joinpath("enrichment.R"),
		cluster_profiler = biorep_dir.joinpath("clusterProfiler.R"),
		A2B_outdir = lambda wildcards: dir2.joinpath(wildcards.spl_pair, wildcards.spl_pair+".A2B"),
		B2A_outdir = lambda wildcards: dir2.joinpath(wildcards.spl_pair, wildcards.spl_pair+".B2A"),
	shell:
		"grep 'A2B' {input.switch_table} > {output.A2B_bed} ;\n"
		
		"bedtools intersect -a {output.A2B_bed} -b {input.gene_bed} -wo | cut -f 10 | sort -u > {output.A2B_gene} ;\n"

		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"    --genelist_file {output.A2B_gene} -o {params.A2B_outdir} -p {wildcards.spl_pair}.A2B ;\n"

		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {output.A2B_gene}"
		"    {params.A2B_outdir}/{wildcards.spl_pair}.A2B ;\n"

		"grep 'B2A' {input.switch_table} > {output.B2A_bed} ;\n"
		
		"bedtools intersect -a {output.B2A_bed} -b {input.gene_bed} -wo | cut -f 10 | sort -u > {output.B2A_gene} ;\n"

		"{params.rscript} {params.cluster_profiler} --ACC2GO {input.acc2go}"
		"    --genelist_file {output.B2A_gene} -o {params.B2A_outdir} -p {wildcards.spl_pair}.B2A ;\n"

		"{params.rscript} {params.enrich} {input.acc2go} {input.kegg2go} {output.B2A_gene}"
		"    {params.B2A_outdir}/{wildcards.spl_pair}.B2A ;\n"

rule switch_plot:
	input:
		chrsize = single_dir.joinpath("01.ref", "genome.chrsize"),
		coolfiles = lambda wildcards: expand(single_dir.joinpath("05.compartment", "01.prepare", 
							"{spl}_{cptmt_res}.cool"), spl=wildcards.spl_pair.split(".minus."), 
							cptmt_res=wildcards.compartment_res),
		gcgenes = lambda wildcards: expand(single_dir.joinpath("05.compartment", "02.PCA", 
							"{spl}_{cptmt_res}", "{spl}_{cptmt_res}.gc_gene.bed"), 
							spl=wildcards.spl_pair.split(".minus."), cptmt_res=wildcards.compartment_res),
	output:
		switch_fig_ = dir2.joinpath("{spl_pair}", "switch_plot.{compartment_res}", "{spl_pair}.chr1.pdf"),
	params:
		python3 = software.get("python3"),
		cptmt_plot = scripts_dir.joinpath("compartment_plot.py"),
		splnames = lambda wildcards: wildcards.spl_pair.split(".minus."),
		outdir = lambda wildcards: dir2.joinpath(wildcards.spl_pair, "switch_plot."+str(wildcards.compartment_res))
	shell:
		"{params.python3} {params.cptmt_plot} --chromsize {input.chrsize} "
		"     --control_name {params.splnames[1]} --control_cool {input.coolfiles[1]}"
		"     --control_cptmt_gcgene {input.gcgenes[1]} --treat_name {params.splnames[0]}"
		"     --treat_cool {input.coolfiles[0]} --treat_cptmt_gcgene {input.gcgenes[0]}"
		"     --pc1_column 4 --outdir {params.outdir} ;\n"

rule cptmt_plot:
	input:
		saddle_dump = expand(single_dir.joinpath("05.compartment", "02.PCA", 
							"{spl}_{compartment_res}", 
							"{spl}_{compartment_res}.saddle.saddledump.npz"), 
							spl=sample_order, compartment_res=compartment_res), 
		cptmt_gcgene = expand(single_dir.joinpath("05.compartment", "02.PCA",
							"{spl}_{compartment_res}", 
							"{spl}_{compartment_res}.gc_gene.bed"), 
							spl=sample_order, compartment_res=compartment_res),
	output:
		saddle_figure = dir2.joinpath("All", "saddle_plot_all.png"),
		cptmt_barfig = dir2.joinpath("All", "compartment_barplot.png"),
	params:
		python3 = software.get("python3"),
		saddle_plot_maker = biorep_dir.joinpath("saddle_plot_maker.py"),
		cptmt_barplot = biorep_dir.joinpath("compartment_barplot.py"),
		sample_names = config.get("sample_order").replace(" ", "").replace(",", " "),
		outpfix = dir2.joinpath("All", "saddle_plot_all"),
		barplot_prefix = dir2.joinpath("All", "compartment_barplot"),
	shell:
		"{params.python3} {params.saddle_plot_maker} --saddledumps {input.saddle_dump}"	
		"    --samplenames {params.sample_names} --outpfix {params.outpfix} ;\n"			
		
		"{params.python3} {params.cptmt_barplot} --names {params.sample_names}"
		"    --compartment {input.cptmt_gcgene} --col 4 --outpfix {params.barplot_prefix} ;\n"

rule diff_cptmt:
	input:
		cptmt_gcgene = lambda wildcards: [single_dir.joinpath("05.compartment",	"02.PCA", 
				f"{spl}_{compartment_res}", f"{spl}_{compartment_res}.gc_gene.bed")
				for spl in wildcards.spl_pair.split(".minus.")],
	output:
		switch_table = dir2.joinpath("{spl_pair}", "{spl_pair}.switch.xls"),
		switch_gc = dir2.joinpath("{spl_pair}", "{spl_pair}.switch_GC.png"),
	params:
		python3 = software.get("python3"),
		rscript = software.get("Rscript"),
		cptmt_gc_gene = biorep_dir.joinpath("cptmt_gc_gene.r"),
		cptmt_switch = biorep_dir.joinpath("cptmt_switch.py"),
		cptmt_switch_plot = biorep_dir.joinpath("cptmt_switch.r"),
		name_pair = lambda wildcards: wildcards.spl_pair.split(".minus."),
		cptmt_outdir = lambda wildcards: dir2.joinpath(wildcards.spl_pair),
		cptmt_dir = lambda wildcards: dir2.joinpath(wildcards.spl_pair),
		outpfix = lambda wildcards: dir2.joinpath(wildcards.spl_pair, f"{wildcards.spl_pair}.switch"),
	shell:
		"{params.python3} {params.cptmt_switch} --control_name {params.name_pair[1]}"
		"    --treat_name {params.name_pair[0]} --control_cptmt {input.cptmt_gcgene[1]}"
		"    --treat_cptmt {input.cptmt_gcgene[0]} --outdir {params.cptmt_outdir} ;\n"
	
		"{params.rscript} {params.cptmt_switch_plot} {output.switch_table} {params.outpfix} ;\n"		

		"{params.rscript} {params.cptmt_gc_gene} --control_name {params.name_pair[1]}"
		"    --treat_name {params.name_pair[0]} --control_cptmt {input.cptmt_gcgene[1]}"
		"    --treat_cptmt {input.cptmt_gcgene[0]} --outdir {params.cptmt_outdir} ;\n"

rule decay_curve:
	input:
		expand(single_dir.joinpath("04.matrix", "{spl}.{genome_res}.cool"), 
							spl=sample_order, genome_res=genome_res),
	output:
		dir1.joinpath("03.diff_pscurve", "decay_curve.png")
	params:
		splnames = sample_order,
		python3 = software.get("python3"),
		contacts_distance = biorep_dir.joinpath("contacts_distance.py"),
		ps_outpfix = dir1.joinpath("03.diff_pscurve", "decay_curve")
	shell:
		"{params.python3} {params.contacts_distance} --in_cool {input}"
		"    --samplename {params.splnames} --outpfix {params.ps_outpfix} ;\n"

rule diff_matrix:
	input: 
		incools = lambda wildcards: [single_dir.joinpath("04.matrix", f"{spl}_{genome_res}",f"{spl}.zScore.matrix.gz") 
						for spl in wildcards.spl_pair.split(".minus.")],
	output:
		genome_figure = dir1.joinpath("01.genome_diffhm", "{spl_pair}", "{spl_pair}.genome.png")
	params:
		python3 = software.get("python3"),
		name_pair = lambda wildcards: wildcards.spl_pair.split(".minus."),
		#diff_matrix = biorep_dir.joinpath("diff_matrix.py"),
		diff_matrix = "/work/frasergen/3D/work/shaojie/script/HiC/hic_prepare/standard_analysis/diff_matrix.py",
		outdir_genome = dir1.joinpath("01.genome_diffhm"),
		outdir_chrom = dir1.joinpath("02.bychr_diffhm"),
		csize=single_dir.joinpath("01.ref","genome.chrsize"),
	run:
		shell("{params.python3} {params.diff_matrix} --matrix genome --wt_mat {input.incools[1]} --case_mat {input.incools[0]} --wt_name {params.name_pair[1]} --case_name {params.name_pair[0]} --outdir {params.outdir_genome} --res {genome_res} -g {params.csize}")
		shell(f"{params.python3} {params.diff_matrix} --matrix chrom --wt_mat {Path(input.incools[1]).parent.joinpath('bychr')} --case_mat {Path(input.incools[0]).parent.joinpath('bychr')} --wt_name {params.name_pair[1]} --case_name {params.name_pair[0]} --outdir {params.outdir_chrom} --res {genome_res} -g {params.csize}")

