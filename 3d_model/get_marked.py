import fire,os
import numpy as np

#pdb文件以N结尾，代表最后一个bin
#开头处有一个N，到N处代表第一个bin
#中间处bin均匀分布

def mark(pdbfile,fai,bed,res,out_f="bin2pdb.txt",out_r="pdb2bin.txt"):
    pdb=np.loadtxt(pdbfile,str)
    Nsite=pdb[np.where(pdb[:,2]=="N")][:,1].astype(int)
    last_bin=[int(i.split()[1])//res+1 for i in open(fai)]
    chains=list(map(chr, (*range(65,91),*range(97,123))))

    b2p,mark=list(),False
    for i in open(bed):
        chrom,start,end=i.split()[:3]
        index,start,end=map(int, (chrom[3:],start,end))
        n1,n2=Nsite[2*index-2:2*index]
        bin_start=start*n1/res if start<res else n2 if start+res>res*last_bin[index-1] else n1+(start-res)/res*(n2-n1-1)/(last_bin[index-1]-2)
        bin_end=end*n1/res if end<res else n2 if end+res>res*last_bin[index-1] else n1+(end-res)/res*(n2-n1-1)/(last_bin[index-1]-2)
        if len(i.split())>3:
            pdb[int(bin_start):int(bin_end),2],mark=i.split()[3],True
        b2p.append("\t".join(map(str,(chrom,start,end,chains[index-1],int(bin_start)+1,int(bin_end)))))
    if mark:
        out_pdb=os.path.basename(pdbfile).replace(".pdb","_marked.pdb")
        np.savetxt(out_pdb,pdb,"%4s%7s%3s%6s%2s%16s%8s%8s%6s%6s")
    with open(out_f,"w") as OUT:
        print(*b2p,sep="\n",file=OUT)

    base_starts=list()
    chroms=list()
    for i in range(pdb.shape[0]):
        index=chains.index(pdb[i,4])
        chroms.append(f"chr{index+1}")
        n,n1,n2=(int(pdb[i,1]),*Nsite[2*index:2*index+2])
        a=(n-1)*res/n1 if n<=n1 else (last_bin[index]-1)*res if n==n2 else (n-n1-1)*(last_bin[index]-2)*res/(n2-n1-1)+res
        base_starts.append(int(a))
    np.savetxt(out_r,np.column_stack((pdb[:,[4,1]],np.array([chroms]).T,np.array([base_starts]).T,pdb[:,5:8])),"%s","\t")

if __name__=="__main__":
    fire.Fire(mark)


