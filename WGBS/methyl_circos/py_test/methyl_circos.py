from pycirclize import Circos
from pycirclize.utils import ColorCycler
import numpy as np
import os,fire
np.random.seed(0)

def get_zscore(data,value):
    data = (data-np.mean(data))/np.std(data)
    return np.clip(data, a_min=-value, a_max=value)

def circos_plot(*txt,fai="genome.fasta.fai",bin_size=100000,outpfix="methyl_test.pdf",z=1):
    gc = np.loadtxt(txt[0],str)
    g_num,cg,chg,chh = map(lambda x:np.loadtxt(x,str)[:,-1].astype(float),txt[1:])
    GC = gc[:,-1].astype(float)

    # init from bed
    fai,bed = [i.split()[:2] for i in open(fai)],fai.replace("fai","bed")
    if not os.path.exists(bed) :
        with open(bed,"w") as OUT:
            print(*("\t".join((i[0],"0",i[1])) for i in fai),sep="\n",file=OUT)
    
    circos = Circos.initialize_from_bed(bed,space=2,start=10,end=350,endspace=False)
    circos.text("Methylation", size=12, r=30) # circos.text位置默认在中心
    
    # get colormap
    ColorCycler.set_cmap("gist_rainbow")
    colors = ColorCycler.get_color_list(len(circos.sectors))
    
    heatmap_cmaps = ("Greens","Oranges","Purples","Reds","Blues")
    texts = ("GC","Genenum","CG","CHG","CHH")
    r_size, r_interval, r_start = 9, 1, 94
    # plot genome and ticks
    unit = int(float(f"1e{len(str(max((i.end for i in circos.sectors))))-1}"))
    for sector, color in zip(circos.sectors, colors):
        sector.text(sector.name.replace("chr",""),size=10,color=color,r=110) # text主要是内容和尺寸，sector.text位置默认在外圈
        genome = sector.add_track((95,100))
        genome.axis(fc=color)
        genome.xticks_by_interval(unit,label_size=5,label_orientation="vertical",                                label_formatter=lambda x:f"{x*1e-6:.0f} Mb")
        genome.xticks_by_interval(unit/2,tick_length=1, show_label=False)
        genome.xticks_by_interval(unit/10,tick_length=0.5, show_label=False)
    
    # 读取bed数据，heatmap plot
        site = np.where(gc[:,0]==sector.name)
        datas= map(lambda x:get_zscore(x[site],z), (GC,g_num,cg,chg,chh))
        
        for idx, cmap, data in zip(range(len(heatmap_cmaps)),heatmap_cmaps,datas):
            r_pos = r_start-(r_size+r_interval)*idx
            r_lim = (r_pos - r_size,r_pos)
            track = sector.add_track(r_lim)
            track.axis(ec="grey")
            track.heatmap(data,width=bin_size,cmap=cmap)
            if sector.name == circos.sectors[0].name:
                circos.text(texts[idx],r=track.r_center,size=8)
                circos.colorbar(bounds=(0.35, 0.6-0.05*idx, 0.3, 0.01),                                       orientation="horizontal", cmap=cmap)

    # 保存
    circos.savefig(outpfix,dpi=750)

if __name__=="__main__":
    fire.Fire(circos_plot)


