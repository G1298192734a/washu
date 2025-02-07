import argparse, os, sys
from os.path import abspath
import numpy as np
from string import Template
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cis_template = """\
<<include etc/housekeeping.conf>>
<<include etc/colors_fonts_patterns.conf>>
data_out_of_range* = trim
<colors>
  <<include etc/brewer.conf>>
</colors>

<ideogram>
  <spacing>
  default = 0.005r
  </spacing>
  radius  = 0.8r
  thickness = 40p
  fill    = yes
  stroke_thickness = 2
  stroke_color = black
  show_bands =yes
  fill_bands =yes
  show_label = yes
  label_font = default
  label_radius = dims(ideogram,radius) + 0.075r
  label_size = 60p
  label_parallel = no
</ideogram>

<image>
  <<include etc/image.conf>>
  dir* = ${outdir}
  file* = ${outfile}
  auto_alpha_steps* = 20
</image>

karyotype = ${karyotype}
chromosomes_display_default = yes
chromosomes_units = 1000000
show_ticks = yes
show_tick_labels=yes

<ticks>
  radius = dims(ideogram,radius_outer)
  multiplier = 1e-6
  <tick>
    spacing = 5u
    size = 8p
    thickness = 2p
    color = black
    show_label = no
    label_size = 20p
    format = %.0f
  </tick>
  <tick>
    spacing = 10u
    size = 16p
    thickness = 2p
    color = black
    show_label = no
    label_size = 40p
    label_offset = 5p
    format = %.0f
  </tick>
</ticks>

<plots>
  z = 0
  <plot>
    type        = histogram
    file        = ${gene_density}
    r1          = 0.98r
    r0          = 0.85r
    min         = 0
    max         = ${max_gene}
    extend_bin  = no
    fill_color  = dred
    color       = dred
    thickness   = 0
  </plot>

  <plot>
    type             = heatmap
    stroke_thickness = 0
    scale_log_base   = 1
    file             = ${link_density}
    r0               = 0.7r
    r1               = 0.78r
    min              = ${link_density_min}
    max              = ${link_density_max}
    color            = ylorrd-9-seq
  </plot>

</plots>

<links>
  z = 5
  radius = 0.69r
  bezier_radius = 0.3r
  bezier_radius_purity = 0.45
  crest = 0.4
  perturb = yes

  <link cis>
    z = 0
    show = yes
    color = vvlgrey
    thickness = 1
    record_limit=25000
    file = ${color_links}
    <rules>
      <rule>
        importance = 205
        condition = var(color) eq 'vvdblue'
        z =90
      </rule>
      <rule>
        importance = 205
        condition = var(color) eq 'vdblue'
        z=50
      </rule>
      <rule>
        importance = 205
        condition = var(color) eq 'dblue'
        z=20
      </rule>
      <rule>
        importance = 205
        condition = var(color) eq 'lblue'
        z=10
      </rule>
    </rules>
  </link>
</links>
<plots>
</plots>"""

def run_circos(myargs):
	outdir = abspath(myargs.outdir)
	try: os.stat(outdir)
	except: os.makedirs(outdir)

	if os.path.getsize(myargs.link_density):
		max_gene = max([int(line.strip().split()[-1]) for line in open(myargs.gene_density)])
		links = np.unique([int(line.strip().split()[-1]) for line in open(myargs.link_density)])
		min_links = int(np.percentile(links, 10))
		max_links = int(np.percentile(links, 90))
		
		cis_info = Template(cis_template)
		cis_conf = cis_info.substitute(karyotype=abspath(myargs.karyotype), outdir=outdir,
					gene_density=abspath(myargs.gene_density), max_gene=max_gene, 
					outfile=f"{myargs.method}_interaction.png",
					link_density=abspath(myargs.link_density), 
					color_links=abspath(myargs.color_links), 
					link_density_min=min_links, link_density_max=max_links)
	
		out_conf = os.path.join(outdir, f"{myargs.method}_interaction.conf")
		with open(out_conf, "w") as cis_buff:
			print(cis_conf, file=cis_buff)
		os.system(f"circos -conf {out_conf}")

	else:
		plt.figure(figsize=(10, 7))
		plt.text(x=0.5, y=0.5, s="No links found", fontsize=18, fontfamily="serif", 
					va="center", ha="center")
		plt.axis('off')
		plt.savefig(os.path.join(outdir, f"{myargs.method}_interaction.png"))
		plt.savefig(os.path.join(outdir, f"{myargs.method}_interaction.svg"))

def get_myargs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-m", "--method", dest="method", required=True, 
		choices=("cis", "trans"), help="generate cis or trans config and run circos [required]")
	parser.add_argument("-k", "--karyotype", dest="karyotype", required=True, 
		help="karyotype file for circos [required]")
	parser.add_argument("-g", "--gene_density", dest="gene_density", required=True, 
		help="gene density file for circos [required]")
	parser.add_argument("-d", "--link_density", dest="link_density", required=True, 
		help="link density file for circos [required]")
	parser.add_argument("-l", "--color_links", dest="color_links", required=True, 
		help="links file contain color info [required]")
	parser.add_argument("-o", "--outdir", dest="outdir", default="./",
		help="output dir path [default:./]")
	return parser.parse_args()

if __name__ == "__main__":
	myargs = get_myargs()
	run_circos(myargs)
