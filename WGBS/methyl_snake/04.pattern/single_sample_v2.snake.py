import yaml, sys, snakemake
from os.path import dirname, abspath
from os.path import join as fjoin
__author__ = "lurui@frasergen.com"
__date__   = "2021.08.18"

snake_dir = "/public/frasergen/3D/pipeline/WGBS/wgbs_pipe/workflow"
configfile: fjoin(dirname(snake_dir), "software.yaml")
configfile: "single_sample.yaml"
script_dir = fjoin(dirname(snake_dir), "scripts")
src_dir = fjoin(dirname(snake_dir), "src")
all_spl = list(config.get('reads').keys())
snakemake.utils.makedirs(all_spl)
work_dir = os.getcwd()

rule all:
    input: 
        expand(abspath(fjoin("{splname}", "{splname}.methylevel.xls")), splname=all_spl),
        expand(abspath(fjoin("{splname}", "{splname}.CG.report.txt.gz")), splname=all_spl),
        expand(abspath(fjoin("{splname}", "{splname}.CG.bw")), splname=all_spl),
        expand(abspath(fjoin("{splname}", "{splname}.mapstat.xls")), splname=all_spl),
        expand(abspath(fjoin("{splname}", "{splname}.genebody.png")), splname=all_spl),
        expand(abspath(fjoin("{splname}", "{splname}.genebody.pdf")), splname=all_spl),
        expand(abspath(fjoin("{splname}", "circos", "{splname}.meth_density.png")), splname=all_spl),
    shell:
        "echo {input}"

rule circos_plot:
    input:
        gene_bed = config.get("gene_bed"),
        cx_report = abspath(fjoin("{splname}", "{splname}.CX_report.txt.gz")),
        reference = abspath(config.get("reference")),
    output:
        methyl_level = abspath(fjoin("{splname}", "circos", "{splname}_meth_level.conf")),
        methyl_density = abspath(fjoin("{splname}", "circos", "{splname}_meth_density.conf")),
        circos_figure = abspath(fjoin("{splname}", "circos", "{splname}.meth_density.png")),
    params:
        mypython = config.get("python3"),
        outdir = abspath(fjoin("{splname}", "circos")),
        methyl_circos = fjoin(script_dir, "methyl_circos.py"),
    shell:
        "[ !-d {params.outdir} ] && mkdir -p {params.outdir} ;\n"
        "{params.mypython} {params.methyl_circos} -cx {input.cx_report}"
        "     --sample {wildcards.splname} --genome {input.reference}"
        "     --genebed {input.gene_bed} --out {params.outdir} --meth_type CG,CHG,CHH ;\n"

        "cd {params.outdir} ;\n"
        "module load circos/0.69-9 ;\n"
        "circos -conf {output.methyl_level} ;\n"
        "circos -conf {output.methyl_density} ;\n"
        "cd {work_dir} ;\n"

rule genebody_enrich:
    input:
        gene_bed = config.get("gene_bed"),
        cg_bw = abspath(fjoin("{splname}", "{splname}.CG.bw")),
        chg_bw = abspath(fjoin("{splname}", "{splname}.CHG.bw")),
        chh_bw = abspath(fjoin("{splname}", "{splname}.CHH.bw")),
    output:
        out_matrix = abspath(fjoin("{splname}", "{splname}.genebody.matrix.gz")),
        out_png = abspath(fjoin("{splname}", "{splname}.genebody.png")),
        out_pdf = abspath(fjoin("{splname}", "{splname}.genebody.pdf")),
    params:
        compute_matrix = fjoin(dirname(config.get("deeptools")), "computeMatrix"),
        plotprofile = fjoin(dirname(config.get("deeptools")), "plotProfile"),
    shell:
        "{params.compute_matrix} scale-regions --numberOfProcessors 16" 
        "     --scoreFileName {input.cg_bw} {input.chg_bw} {input.chh_bw}"
        "     --regionsFileName {input.gene_bed} --beforeRegionStartLength 3000"
        "     --afterRegionStartLength 3000 --regionBodyLength 3000"
        "     -o {output.out_matrix} --samplesLabel CG CHG CHH ;\n"

        """{params.plotprofile} --matrixFile {output.out_matrix} --dpi 600 \
             --outFileName {output.out_png} --perGroup --plotWidth 20 \
             --plotHeight 13 --colors '#0fbe3c' '#be3c0f' '#3c0fbe' ;\n """

        """ {params.plotprofile} --matrixFile {output.out_matrix} --dpi 600 \
             --outFileName {output.out_pdf} --perGroup --plotWidth 20 \
             --plotHeight 13 --colors '#0fbe3c' '#be3c0f' '#3c0fbe' ;\n """

rule output_stat:
    input:
        chromsize = abspath(config.get("chromsize")),
        align_stat = abspath(fjoin("{splname}", "{splname}_PE_report.txt")),
        mbias_stat = abspath(fjoin("{splname}", "{splname}.deduplicated.M-bias.txt")),
        dedup_stat = abspath(fjoin("{splname}", "{splname}_pe.deduplication_report.txt")),
        dedup_report = abspath(fjoin("{splname}", "{splname}.CX_report.txt.gz")),
    output:
        map_stat = abspath(fjoin("{splname}", "{splname}.mapstat.xls")),
        mbias_pdf = abspath(fjoin("{splname}", "{splname}.mbias.pdf")),
        cover_stat = abspath(fjoin("{splname}", "{splname}.coverage_stat.xls")),
        methyl_stat = abspath(fjoin("{splname}", "{splname}.methylevel.xls")),
        cbase_depth = abspath(fjoin("{splname}", "{splname}.cbase_depth.pdf")),
        cpoint_mrate = abspath(fjoin("{splname}", "{splname}.cpoint_mrate.pdf")),
        cg_bw = abspath(fjoin("{splname}", "{splname}.CG.bw")),
        chg_bw = abspath(fjoin("{splname}", "{splname}.CHG.bw")),
        chh_bw = abspath(fjoin("{splname}", "{splname}.CHH.bw")),
    params:
        bedg2bw = fjoin(src_dir, "bedGraphToBigWig"),
        bedtools = config.get("bedtools"),
        outpfix = abspath(fjoin("{splname}", "{splname}")),
        rscript = config.get("rscript"),
        mypython = config.get("python3"),
        bismark_stat = fjoin(script_dir, "bismark_stat.r"),
        methyl_ratio = fjoin(script_dir, "methyl_ratio.py"),
    shell:
        "{params.rscript} {params.bismark_stat} --mbias_stat {input.mbias_stat}"
        "     --align_stat {input.align_stat} --dedup_stat {input.dedup_stat}"
        "     --samplename {wildcards.splname} --outdir {wildcards.splname};\n"
    	
        "{params.mypython} {params.methyl_ratio} --cx_report {input.dedup_report}"
        "     --chromsize {input.chromsize} --bedg2bw {params.bedg2bw}"
        "     --bedtools {params.bedtools} --outpfix {params.outpfix} ;\n"

rule split_cx:
    input:
        dedup_report = abspath(fjoin("{splname}", "{splname}.CX_report.txt.gz")),
    output:
        out_report = abspath(fjoin("{splname}", "{splname}.CG.report.txt.gz")),
    params:
        outdir = abspath("{splname}"),
        mypython = config.get("python3"),
        methyl_utils = fjoin(script_dir, "methyl_utils.py"),
    shell:
        "{params.mypython} {params.methyl_utils} split_cx"
        "    --infile {input.dedup_report} --outdir {params.outdir} ;\n"

rule methyl_report:
    input:
        align_stat = abspath(fjoin("{splname}", "{splname}_PE_report.txt")),
        nucl_stat = abspath(fjoin("{splname}", "{splname}_pe.nucleotide_stats.txt")),
        dedup_bam = abspath(fjoin("{splname}", "{splname}.deduplicated.bam")),
        dedup_stat = abspath(fjoin("{splname}", "{splname}_pe.deduplication_report.txt")),
    output:
        cg_dedup = abspath(fjoin("{splname}", "CpG_context_{splname}.deduplicated.txt.gz")),
        chg_dedup = abspath(fjoin("{splname}", "CHG_context_{splname}.deduplicated.txt.gz")),
        chh_dedup = abspath(fjoin("{splname}", "CHH_context_{splname}.deduplicated.txt.gz")),
        mbias_stat = abspath(fjoin("{splname}", "{splname}.deduplicated.M-bias.txt")),
        split_stat = abspath(fjoin("{splname}","{splname}.deduplicated_splitting_report.txt")),
        dedup_bdg = abspath(fjoin("{splname}", "{splname}.deduplicated.bedGraph.gz")),
        dedup_cov = abspath(fjoin("{splname}", "{splname}.deduplicated.bismark.cov.gz")),
        dedup_report = abspath(fjoin("{splname}", "{splname}.CX_report.txt.gz")),
        final_stat = abspath(fjoin("{splname}", "{splname}.integrate.html")),
    params:
        samtools = config.get("samtools"),
        outbedg = "{splname}.deduplicated.bedGraph.gz",
        out_cstat = "{splname}.dedup.CX_report.txt.gz",
        ref_dir = dirname(config.get("reference")),
        cov_to_c = fjoin(dirname(config.get("bismark")), "coverage2cytosine"),
        to_bedg = fjoin(dirname(config.get("bismark")), "bismark2bedGraph"),
        final_result = fjoin(dirname(config.get("bismark")), "bismark2report"),
        extractor = fjoin(dirname(config.get("bismark")), "bismark_methylation_extractor"),
    threads: 32
    shell:
        "{params.extractor} --paired-end --no_overlap --comprehensive --gzip"
        "     --parallel {threads} --buffer_size 120G --ignore 5 --ignore_r2 5"
        "     --output {wildcards.splname} {input.dedup_bam} ;\n"
    
        "{params.to_bedg} --CX --buffer_size 120G --dir {wildcards.splname}"
        "     --output {params.outbedg} {output.cg_dedup} {output.chg_dedup} {output.chh_dedup} ;\n"

        "{params.cov_to_c} --CX --gzip --genome_folder {params.ref_dir} --output {wildcards.splname}"
        "     --dir {wildcards.splname} {output.dedup_cov} ;\n" 

        "{params.final_result} --alignment_report {input.align_stat} --dedup_report {input.dedup_stat}"
        "     --splitting_report {output.split_stat} --mbias_report {output.mbias_stat}"
        "     --nucleotide_report {input.nucl_stat} -o {output.final_stat} ;\n" 

rule alignment:
    input:
        reference = abspath(config.get("reference")),
        r1 = lambda wildcards: config["reads"][wildcards.splname][0],
        r2 = lambda wildcards: config["reads"][wildcards.splname][1],
    output:
        raw_bam = abspath(fjoin("{splname}", "{splname}_pe.bam")),
        nucl_stat = abspath(fjoin("{splname}", "{splname}_pe.nucleotide_stats.txt")),
        align_stat = abspath(fjoin("{splname}", "{splname}_PE_report.txt")),
        dedup_bam = abspath(fjoin("{splname}", "{splname}.deduplicated.bam")),
        dedup_stat = abspath(fjoin("{splname}", "{splname}_pe.deduplication_report.txt")),
    params:
        bismark = config.get("bismark"),
        bowtie2_dir = dirname(config.get("bowtie2")),
        samtools = dirname(config.get("samtools")),
        ref_dir = dirname(config.get("reference")),
        deduplicate = fjoin(dirname(config.get("bismark")), "deduplicate_bismark"),
    threads: 32
    shell:
        "{params.bismark} --nucleotide_coverage --path_to_bowtie2 {params.bowtie2_dir}"
        "     --output_dir {wildcards.splname} --gzip --bam --samtools_path {params.samtools}"
        "     --basename {wildcards.splname} -p {threads} --temp_dir {wildcards.splname}"
        "     {params.ref_dir} -1 {input.r1} -2 {input.r2} ;\n"
        
        "{params.deduplicate} --paired --bam --samtools_path {params.samtools}"
        "     --output_dir {wildcards.splname} --outfile {wildcards.splname} {output.raw_bam} ;\n"
