from pandas import read_table
import seaborn as sns
import matplotlib.pyplot as plt
import fire

def plot_coverage(depth):
  plot_data=read_table(depth,dtype={"0":"str","4":"float","5":"float"},usecols=[0,4,5],header=None,engine="c",na_filter=False,low_memory=False,comment="#").iloc[:-1]

  fig, axs = plt.subplots(1, 1)
  sns.barplot(x=plot_data[0],y=plot_data[5],ax=axs)
  axs.set_xlabel("")
  axs.set_ylabel("Mean Depth")
  axs.set_title("Chromosome Coverage and Depth")

  axs1=axs.twinx()
  sns.scatterplot(x=plot_data[0],y=plot_data[4],marker="o",edgecolor="r",ax=axs1,facecolor="none",s=20,linewidth=1)
  sns.lineplot(x=plot_data[0],y=plot_data[4],color="r",ls="-.",ax=axs1,linewidth=0.5)
  axs1.set_ylabel("Prop of Coverage")
  axs1.set_ylim(0,100)
  plt.setp(axs.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
  plt.savefig("coverage_distribution_chr.pdf")

if __name__=="__main__":
  fire.Fire(plot_coverage)


