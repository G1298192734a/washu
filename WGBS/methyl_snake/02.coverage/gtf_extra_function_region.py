import pybedtools as pb
import fire,re

def gtf_extract(gtf,fai):
  gtf = pb.BedTool(gtf).saveas()
  fai = {i.split()[0]:(0,int(i.split()[1])) for i in open(fai)}
  featuretypes = ("transcript","exon","CDS")
  def subset_featuretypes(featuretype):
    return gtf.filter(
             lambda x,y: True if x[2]==y else False, featuretype
           ).sort().saveas()
  rnas, exons, cds = map(subset_featuretypes, featuretypes)

  pos_gene=rnas.filter(
                 lambda x: True if x[6]=="+" else False
               ).merge().saveas()
  neg_gene=rnas.filter(
                 lambda x: True if x[6]=="-" else False
               ).merge().saveas()
  pos_promoter=pos_gene.slop(l=2000,r=0,g=fai).saveas()
  pos_downstream=pos_gene.slop(l=0,r=2000,g=fai).saveas()
  neg_promoter=neg_gene.slop(l=0,r=2000,g=fai).saveas()
  neg_downstream=neg_gene.slop(l=2000,r=0,g=fai).saveas()
  promoters = pos_promoter.cat(neg_promoter)\
                .subtract(rnas).remove_invalid().saveas()
  downstreams = pos_downstream.cat(neg_downstream)\
                .subtract(rnas).remove_invalid().saveas()
  introns = rnas.merge().subtract(exons).remove_invalid()

  utrs = exons.merge().subtract(cds).remove_invalid().saveas()
  inter=utrs.intersect(rnas,wb=True)
  fifutrs=inter.filter(
            lambda x:True if (x[9]=="+" and int(x[1])-int(x[6])<=int(x[7])-int(x[2])) or (x[9]=="-" and int(x[7])-int(x[2])<=int(x[1])-int(x[6])) else False
          ).each(lambda x: x[:3])
  triutrs=inter.filter(
            lambda x:True if (x[9]=="+" and int(x[1])-int(x[6])>int(x[7])-int(x[2])) or (x[9]=="-" and int(x[7])-int(x[2])>int(x[1])-int(x[6])) else False
          ).each(lambda x: x[:3])

  with open("genebody_region.bed","w") as OUT:
    for i,j in zip(("upstream","genebody","downstream"),(promoters,rnas.merge(),downstreams)):
      print(*("\t".join((*q,i+str(p),".",".")) for p,q in enumerate(j)),sep="\n",file=OUT)
  with open("function_region.bed","w") as OUT:
    for i,j in zip(("promoter","5utr","exon","intron","3utr","downstream"),(promoters,fifutrs,exons.merge(),introns,triutrs,downstreams)):
      print(*("\t".join((*q,i+str(p),".",".")) for p,q in enumerate(j)),sep="\n",file=OUT)

  #for i,j in zip(("exon","intron","utr"),(exons,introns,utrs)):
    #with open(f"{i}.bed","w") as OUT:
      #print(*("\t".join((*q,i+str(p),".",".")) for p,q in enumerate(j)),sep="\n",file=OUT)
  pb.cleanup(verbose=False)

if __name__=="__main__":
  fire.Fire(gtf_extract)


