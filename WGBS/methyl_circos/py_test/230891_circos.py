from pycirclize import Circos
from pycirclize.utils import ColorCycler
from pycirclize.parser import Gff
import numpy as np
np.random.seed(0)

# Initialize Circos from BED chromosomes
circos = Circos.initialize_from_bed("carguta.fasta.bed", space=2, start=10, end=350, endspace=False)
gff = Gff("/work/frasergen/backup/3d/project/Interactome/230891_ctinidia_arguta_20240204/00.ref/carguta.gff")

red_fpkm = np.loadtxt("/work/frasergen/backup/3d/project/Interactome/230891_ctinidia_arguta_20240204/02.RNA_result/01.genome_gene_express/RED.fpkm4washu",str)
red_fpkm = red_fpkm[np.where(red_fpkm[:,3]!="0")][:,:4]
r_pos,r_fpkm=np.mean(red_fpkm[:,1:3].astype(float),axis=1),np.log(red_fpkm[:,-1].astype(float))
r_vmax,r_p_fpkm,r_n_fpkm = np.max(r_fpkm),np.where(r_fpkm>=0,r_fpkm,0),np.where(r_fpkm<0,r_fpkm,0)

green_fpkm = np.loadtxt("/work/frasergen/backup/3d/project/Interactome/230891_ctinidia_arguta_20240204/02.RNA_result/01.genome_gene_express/GREEN.fpkm4washu",str)
green_fpkm = green_fpkm[np.where(green_fpkm[:,3]!="0")][:,:4]
g_pos,g_fpkm=np.mean(green_fpkm[:,1:3].astype(float),axis=1),np.log(green_fpkm[:,-1].astype(float))
g_vmax,g_p_fpkm,g_n_fpkm = np.max(g_fpkm),np.where(g_fpkm>=0,g_fpkm,0),np.where(g_fpkm<0,g_fpkm,0)

red_compt = np.loadtxt("/work/frasergen/backup/3d/project/Interactome/230891_ctinidia_arguta_20240204/03.HiC/05.compartment/02.PCA/RED_100000/RED_100000.compartment.bed",str)
green_compt = np.loadtxt("/work/frasergen/backup/3d/project/Interactome/230891_ctinidia_arguta_20240204/03.HiC/05.compartment/02.PCA/GREEN_100000/GREEN_100000.compartment.bed",str)

red_tad = np.loadtxt("/work/frasergen/backup/3d/project/Interactome/230891_ctinidia_arguta_20240204/03.HiC/06.tad/02.cooltools/RED.40000/TAD_genome.bed",str)
green_tad = np.loadtxt("/work/frasergen/backup/3d/project/Interactome/230891_ctinidia_arguta_20240204/03.HiC/06.tad/02.cooltools/GREEN.40000/TAD_genome.bed",str)

link = np.loadtxt("/work/frasergen/backup/3d/project/Interactome/230891_ctinidia_arguta_20240204/06.personal/analysis20240605_circos/link_color.new.circos",str)[:,:6]

# Create chromosome color mapping
ColorCycler.set_cmap("gist_rainbow")
chr_names = [i.name for i in circos.sectors]
colors = ColorCycler.get_color_list(len(chr_names))
chr_name2color = {name: color for name, color in zip(chr_names, colors)}
gene_feats = gff.get_seqid2features("gene")

# Plot chromosome name & xticks
unit = int(float(f"1e{len(str(max((i.end for i in circos.sectors))))-1}"))
for sector in circos.sectors:
  sector.text(sector.name, r=120, size=10, color=chr_name2color[sector.name])
  genome = sector.add_track((96,100))
  genome.axis()
  if sector.name == circos.sectors[0].name:
    circos.text("ctinidia arguta",r=genome.r_center+4,size=10)
#    circos.line(r=genome.r_center,color="red",lw=3)
#    circos.text("",r=genome.r_center,size=4)
  genome.xticks_by_interval(unit,label_size=8,label_orientation="vertical",                                 label_formatter=lambda x:f"{x*1e-6:.0f} Mb")
  # Plot forward & reverse gene
  for feat in gene_feats.get(sector.name):
    if feat.location.strand == 1:
      genome.genomic_features([feat], r_lim=(98, 100), fc="red")
    else:
      genome.genomic_features([feat], r_lim=(96, 98), fc="blue")

  # RNA fpkm
  rna = sector.add_track((85,96))
  rna.axis()
  if sector.name == circos.sectors[0].name:
    circos.text("RNA FPKM",r=rna.r_center,size=8)
  site=np.where(red_fpkm[:,0]==sector.name)
  rna.fill_between(r_pos[site], r_p_fpkm[site], 0, vmin=-r_vmax, vmax=r_vmax, color="red")
  rna.fill_between(r_pos[site], r_n_fpkm[site], 0, vmin=-r_vmax, vmax=r_vmax, color="grey")
  rna = sector.add_track((76,85))
  rna.axis()
  site=np.where(green_fpkm[:,0]==sector.name)
  rna.fill_between(g_pos[site], g_p_fpkm[site], 0, vmin=-g_vmax, vmax=g_vmax, color="green")
  rna.fill_between(g_pos[site], g_n_fpkm[site], 0, vmin=-g_vmax, vmax=g_vmax, color="grey")

  # compartment
  compt = sector.add_track((72,76))
  compt.axis()
  if sector.name == circos.sectors[0].name:
    circos.text("Compartment",r=compt.r_center,size=8)
  for i in red_compt[np.where(red_compt[:,0]==sector.name)]:
    if i[-1] == "A":
      compt.rect(int(i[1]),int(i[2]),r_lim=(74,76), fc="red")
    else:
      compt.rect(int(i[1]),int(i[2]),r_lim=(72,74), fc="blue")
  compt = sector.add_track((68,72))
  compt.axis()
  for i in green_compt[np.where(green_compt[:,0]==sector.name)]:
    if i[-1] == "A":
      compt.rect(int(i[1]),int(i[2]),r_lim=(68,70), fc="red")
    else:
      compt.rect(int(i[1]),int(i[2]),r_lim=(70,72), fc="blue")
    
  # tad
  tad = sector.add_track((65,68))
  tad.axis()
  if sector.name == circos.sectors[0].name:
    circos.text("TAD",r=tad.r_center,size=8)
  data,r_lims=red_tad[np.where(red_tad[:,0]==sector.name)],((67,68),(66,67),(65,66))
  for i,j in zip(data,(0,1,2,1)*(data.shape[0]//4+1)):
    tad.rect(int(i[1]),int(i[2]),r_lim=r_lims[j], fc="red")
  tad = sector.add_track((62,65))
  tad.axis()
  data,r_lims=green_tad[np.where(green_tad[:,0]==sector.name)],((62,63),(63,64),(64,65))
  for i,j in zip(data,(0,1,2,1)*(data.shape[0]//4+1)):
    tad.rect(int(i[1]),int(i[2]),r_lim=r_lims[j], fc="green")

  # link
  for i in link:
    reg1,reg2 = (i[0].lower(),int(i[1]),int(i[2])),(i[3].lower(),int(i[4]),int(i[5]))
    if i[0]=="Chr1":
      circos.link(reg1,reg2,r1=61,r2=61,color="#5F9EA0")

# save
circos.savefig("circos.png", dpi=300)

