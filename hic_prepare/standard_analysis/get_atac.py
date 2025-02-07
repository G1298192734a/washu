import os,fire
from pathlib import Path
import shutil

def atac_pre(mapfile,atac):
    atac2hic=dict((i.split()[::-1] for i in open(mapfile).readlines()[1:]))
    ana=Path(atac).resolve()
    
    shutil.copy(next(ana.glob("02*/00*/merge_peak.fa")),"merge_peak.fa")
    for i,j in atac2hic.items():
        if "-vs-" in i: continue
        if not os.path.exists(os.path.join("bw",j)):
            os.mkdir(os.path.join("bw",j))
        for k in ana.glob(f"01*/01*/*/{i}*CPM.bw"):
            os.symlink(k,os.path.join("bw",j,k.name.replace("CPM.","")))

    for i,j in atac2hic.items():
        if "-vs-" not in i: continue
        a,b,c,d=(*i.split("-vs-"),*j.split(".minus."))

        merge=next(ana.glob(f"02*/00*/group_peak/{a}*bed"))
        shutil.copy(merge,merge.name.replace(a,d))
        merge=next(ana.glob(f"02*/00*/group_peak/{b}*bed"))
        shutil.copy(merge,merge.name.replace(b,c))

        deg=next(ana.glob(f"02*/01*/{i}/{i}.D*tsv"))
        shutil.copy(deg,f"{j}_diff.tsv")

if __name__=="__main__":
    fire.Fire(atac_pre)

