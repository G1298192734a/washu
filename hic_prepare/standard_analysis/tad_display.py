import fanc, fire, re
import pysam,pyBigWig
from pathlib import Path
from bisect import bisect_right
import numpy as np
import seaborn as sns
import pybedtools as pb
import fanc.plotting as kplot
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
plt.rcParams["font.size"]=18

def tad_plot(cool,bw,bed,region,outpref,hmax):
	bigwig=pyBigWig.open(bw)
	cool,bw,bed,region = map(fanc.load,(cool,bw,bed,region))
	for i in region:
		overlap=list(bed.intersect(pb.BedTool("\t".join(i),from_string=True),wa=True))
		target=f"{i.chrom}:{overlap[0].start}-{overlap[-1].end}"
		if overlap[-1].end-overlap[0].start>10000000:continue
		m_dist=overlap[-1].end-overlap[0].start
		mj=10**len(str(m_dist))
		mj,mn=mj//5,mj//10

		tad_cool=kplot.TriangularMatrixPlot(cool,vmax=hmax,max_dist=m_dist//2,colormap=sns.blend_palette(["white", "red"], as_cmap=True),draw_x_axis=False,ylabel=i.chrom,aspect=0.3)
		vmax,vmin=map(lambda x:bigwig.stats(i.chrom,overlap[0].start,overlap[-1].end,type=x)[0],("max","min"))
		tad_bw=kplot.LinePlot(bw,style='mid',fill=False,colors="green",plot_kwargs={'alpha':1},show_legend=False,axes_style="ticks",aspect=0.1,ylim=(-1 if vmin=="Null" else vmin,1 if vmax=="Null" else vmax),draw_minor_ticks=True)
		gfig=kplot.GenomicFigure([tad_cool,tad_bw],ticks_last=True,width=8)
		fig,axes=gfig.plot(target)
		axes[1].xaxis.set_major_locator(plt.MultipleLocator(mj))
		axes[1].xaxis.set_minor_locator(plt.MultipleLocator(mn))
		axes[1].axline((0,0), slope=0, color="grey", alpha=0.5, lw=1, ls="--")
		for j in overlap:
			high=(j.end-j.start)//2
			axes[0].add_patch(Polygon(xy=np.array([[j.start,0],[j.start+high,high],[j.end,0]]),closed=False,lw=1.5,ec="black",fc="none"))
		fig.savefig(f"{outpref}.{i.chrom}_{overlap[0].start}_{overlap[-1].end}.pdf", dpi=750)
		fig.savefig(f"{outpref}.{i.chrom}_{overlap[0].start}_{overlap[-1].end}.png", dpi=750)

def main(inref, binsize, cool_matrix, is_bw, tad_bed, outdir, point_num=10, vmax=0.01):
	Path(outdir).mkdir(parents=True, exist_ok=True)
	
	infa = pysam.FastaFile(inref)
	falen = dict(zip(infa.references, infa.lengths))
	falen = sorted(falen.items(), key=lambda x:x[1], reverse=True)[: point_num]
	falen = dict(falen)

	recomp = re.compile(r"[ATCG]+")
	outfile = Path(outdir).joinpath("plot_region.bed")

	with open(outfile, "w") as outbuff:
		for fa in pysam.FastxFile(inref):
			if fa.name not in falen: continue
			seq = fa.sequence.upper()
			start, end = 0, 0
	
			for record in re.finditer(recomp, seq):
				gap_start = min(record.start(), record.end())
				gap_end = max(record.start(), record.end())
				if gap_end - gap_start > end - start:
					start, end = gap_start, gap_end

			mid = (start + end) // 2 
			plot_start = max(0, mid - binsize * 100)
			plot_end = min(mid + binsize * 100, falen[fa.name])
			print(fa.name, plot_start, plot_end, sep="\t", file=outbuff)
	spl_res = Path(outdir).parent.name
	tad_plot(cool_matrix, is_bw, tad_bed, outfile, Path(outdir).joinpath(spl_res), vmax)

if __name__ == '__main__':
	fire.Fire(main)


