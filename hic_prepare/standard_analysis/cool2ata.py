#!/public/frasergen/PUB/software/python/Python-3.9.12/bin/python3
from coolpuppy import coolpup
import cooler,fire
import pandas as pd

def filter_bed(df,minsize=60000,maxsize=3000000):
  length=df["end"]-df["start"]
  return df[(length >= minsize) & (length <= maxsize)]

def cool2ata(Type,cool,bed,outata):
  func={
    #"domain": lambda clr,sites: coolpup.pileup(clr,filter_bed(sites),features_format="bed",nshifts=10,local=True,rescale=True,nproc=32,seed=100,min_diag=0),
    "domain": lambda clr,sites: coolpup.pileup(clr,sites,features_format="bed",nshifts=10,local=True,rescale=True,nproc=32,seed=100,min_diag=0),
    "border": lambda clr,sites: coolpup.pileup(clr,sites,features_format="bed",nshifts=10,local=True,min_diag=0,seed=100,flank=500000),
    "loop"  : lambda clr,sites: coolpup.pileup(clr,sites,features_format="bedpe",nshifts=10),
    "tad"   : lambda clr,sites: coolpup.pileup(clr,sites,features_format="bed",nshifts=10,local=True,rescale=True,seed=100,min_diag=0),
  }

  clr=cooler.Cooler(cool)
  sites=pd.read_table(bed,header=None)
  sites.columns=("chrom","start","end") if sites.shape[1]==3 else ("chrom1","start1","end1","chrom2","start2","end2")

  ata=func[Type](clr,sites)
  comment=f"# baselist: {bed}"
  with open(outata,"w") as OUT:
    print(comment,file=OUT)
    pd.DataFrame(ata.loc[0, "data"]).dropna().to_csv(OUT,sep="\t",header=False,index=False)

if __name__=="__main__":
  fire.Fire(cool2ata)


