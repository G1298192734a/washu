import os,fire
from pathlib import Path
import shutil

def rna_pre(mapfile,rna,bam,fpkm="01.genome_gene_express",deg="02.diff_gene_express"):
    hic2rna=dict((i.split() for i in open(mapfile).readlines()[1:]))
    rna=Path(rna).resolve()

    matrix=list()
    with open(next(rna.glob("02*/gene.fpkm_matrix*"))) as OUT:
        a=OUT.readline().split()
        start,end=map(a.index, ("gene_id","Gene_name"))
        matrix=[a,*(i.split()[start:end] for i in OUT)]

    for i,j in hic2rna.items():
        if ".minus." in i: continue
        if not os.path.exists(os.path.join("03.BAM_2_coverage",i)):
            os.mkdir(os.path.join("03.BAM_2_coverage",i))
        for k in Path(bam).rglob(f"{j}*bam*"):
            a=os.path.join("03.BAM_2_coverage",i,f"{k.name.split('.')[0]}.bam")
            if os.path.exists(a): continue
            os.symlink(k,a)

        sites=[matrix[0].index(k) for k in matrix[0] if k.startswith(j)]
        f=lambda x: "\t".join((x[k] for k in (0,*sites)))
        with open(os.path.join(fpkm,f"{i}_gene_fpkm.xls"),"w") as OUT:
            print(*map(f,matrix),sep="\n",file=OUT)

    for i,j in hic2rna.items():
        if ".minus." not in i: continue
        a=os.path.join(deg,f"{i}.known.DEG.xls")
        shutil.copy(next(rna.glob(f"06*/genes/*known.DEG*")),a)
        a=os.path.join(deg,f"{i}.DEseq2.before-filter.xls")
        shutil.copy(next(rna.glob(f"*/genes/{j}*before-filter*")),a)
        
        reps=list()
        for k in i.split(".minus."):
            f=lambda x:x.startswith(hic2rna[k])
            reps+=(*map(lambda x:f"{k}\t{x}", filter(f,matrix[0])),)
        with open(os.path.join(deg,f"{i}_biorep.txt"),"w") as OUT:
            print(*reps,sep="\n",file=OUT)

if __name__=="__main__":
    fire.Fire(rna_pre)


