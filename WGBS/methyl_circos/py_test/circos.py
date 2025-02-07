from pycirclize import Circos
from pycirclize.utils import ColorCycler
import matplotlib.pyplot as plt
import numpy as np
import fire

def get_median(data):
    data = (data-np.mean(data))/np.std(data)
    p25,p75=np.percentile(data,25),np.percentile(data,75)
    data[np.where(data>p75)],data[np.where(data<p25)]=p75,p25
    return data

def circos_plot(fai,outpfix,*txt,bin_size = 100000):
  gc = np.loadtxt(txt[0],str)
  g_num,cg,chg,chh = map(get_median,map(lambda x:np.loadtxt(x,str)[:,-1].astype(float),txt[1:]))
  
  fai=np.loadtxt("genome.fasta.fai",str)[:,:2]
  sectors = {i:j for i,j in zip(fai[:,0],fai[:,1].astype(int)//bin_size)}
  unit = int(float(f"1e{len(str(max(sectors.values())))-1}"))
  circos = Circos(sectors, space=2, start=15,endspace=False)
  ColorCycler.set_cmap("gist_rainbow")
  colors = ColorCycler.get_color_list(len(circos.sectors))
  for sector,color in zip(circos.sectors,colors):
    site = np.where(gc[:,0]==sector.name)
    # Plot sector name
    sector.text(sector.name, r=110, size=15)
  
    # Plot xtick
    x_track = sector.add_track((95, 100), r_pad_ratio=0.1)
    x_track.axis(fc=color)
    x_track.xticks_by_interval(unit,label_orientation="vertical")
    x_track.xticks_by_interval(unit/2, tick_length=1, show_label=False)
    x_track.xticks_by_interval(unit/10, tick_length=0.5, show_label=False)
    if sector==circos.sectors[-1]:x_track.yticks([0.5], ["Unit: Mb"], vmin=0, vmax=1)
  
    # Plot GC
    track1 = sector.add_track((85, 94))
    track1.axis()
    track1.heatmap(get_median(gc[site][:,-1].astype(float)),cmap="Greens")
    if sector==circos.sectors[-1]:track1.yticks([0.5], ["GC"], vmin=0, vmax=1)
  
    # Plot genenum
    track2 = sector.add_track((75, 84))
    track2.axis()
    track2.heatmap(g_num[site],cmap="Oranges")
    if sector==circos.sectors[-1]:track2.yticks([0.5], ["Genenum"], vmin=0, vmax=1)
  
    # Plot CG,CHG,CHH
    track3 = sector.add_track((65, 74))
    track3.axis()
    track3.heatmap(cg[site],cmap="Purples")
    if sector==circos.sectors[-1]:track3.yticks([0.5], ["CG"], vmin=0, vmax=1)
    track4 = sector.add_track((55, 64))
    track4.axis()
    track4.heatmap(chg[site],cmap="Reds")
    if sector==circos.sectors[-1]:track4.yticks([0.5], ["CHG"], vmin=0, vmax=1)
    track5 = sector.add_track((45, 54))
    track5.axis()
    track5.heatmap(chh[site],cmap="Blues")
    if sector==circos.sectors[-1]:track5.yticks([0.5], ["CHH"], vmin=0, vmax=1)
  
  # Plot text description
#  text_common_kws = dict(ha="left", va="center", size=8)
#  circos.text("Unit: Mb", r=97.5, color="black", **text_common_kws)
#  circos.text("GC",    r=89, color="black", **text_common_kws)
#  circos.text("Genenum", r=79, color="black", **text_common_kws)
#  circos.text("CG", r=69, color="black", **text_common_kws)
#  circos.text("CHG", r=59, color="black", **text_common_kws)
#  circos.text("CHH", r=49, color="black", **text_common_kws)
  
  # Plot colorbar
#  circos.colorbar(bounds=(0.35, 0.55, 0.3, 0.01), orientation="horizontal")
#  circos.colorbar(bounds=(0.35, 0.45, 0.3, 0.01), orientation="horizontal", cmap="viridis")

  circos.plotfig()
  plt.savefig(f"{outpfix}.pdf", format='pdf', dpi=750, bbox_inches='tight')

if __name__=="__main__":
    fire.Fire(circos_plot)


