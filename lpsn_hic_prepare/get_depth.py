import fire,sys

def get_depth(size,rep,out):
    a=[1e-6 for i in range(size)]
    for i in sys.stdin:
        site,value=map(int,i.split()[1:])
        a[site-1]+=value
    with open(out,"w") as OUT:
        print(*(f"chr1\t{i+1}\t{a[i]/int(rep)}" for i in range(size)),sep="\n",file=OUT)

if __name__=="__main__":
    fire.Fire(get_depth)

