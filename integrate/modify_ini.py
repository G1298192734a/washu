from configparser import ConfigParser
from collections import OrderedDict,defaultdict
from pybedtools import BedTool
import pyBigWig as pb
import numpy as np
import fire,re

class MultiDict(OrderedDict):
  """
  Class to allow identically named
  sections in configuration file
  by appending the section number as
  for example:
  1. section name
  """
  _unique = 0

  def __setitem__(self, key, val):
    if isinstance(val, OrderedDict):
      self._unique += 1
      key = f"{str(self._unique)}.{key}"
    OrderedDict.__setitem__(self, key, val)

def modify(ini,region):
  bed=BedTool(region)
  site=[bed[0].chrom,bed[0].start,bed[0].end]
  parser = ConfigParser(dict_type=MultiDict, strict=False)
  parser.read_file(open(ini, 'r'))

  bw_dict,insu_sec,insu_max=defaultdict(list),list(),list()
  for i in parser.sections():
    if i.endswith("bw"): 
      bw_max=pb.open(parser.get(i,"file")).stats(*site,exact=True,type="max")
      bw_dict[parser.get(i,"title").split()[-1]].append((i,*bw_max))
    elif re.search("insulation",i):
      insu_sec.append(i)
      insu_max+=pb.open(parser.get(i,"file")).stats(*site,exact=True,type="max")
      insu_max+=pb.open(parser.get(i,"file")).stats(*site,exact=True,type="min")
    elif re.search("gff genes",i):
      genes = BedTool(parser.get(i,"file")).filter(lambda x:re.search("^(transcript|mRNA)$",x[2]))
      label_num=len(genes.intersect(bed))
      parser.set(i,"max_labels",str(label_num))
      parser.set(i,"height",str(0.5 if label_num==1 else np.log2(label_num)))

  for i,j in bw_dict.items():
    vmax=str(max(p[1] for p in j))
    for p in j:
      parser.set(p[0],"max_value",vmax)
  insu_max=[j for j in insu_max if isinstance(j,float)]
  if insu_max:
    parser.set(insu_sec[0],"max_value",str(max(insu_max)))
    parser.set(insu_sec[0],"min_value",str(min(insu_max)))

  with open(re.sub("ini$","modi.ini",ini),"w") as output:
    for i in parser.sections():
      print("[{0}]".format(re.sub("^\d+.","",i)),file=output)
      print(*(re.sub("\n","\n  "," = ".join(j)) for j in parser.items(i)),sep="\n",file=output)

if __name__=="__main__":
  fire.Fire(modify)


