import numpy as np
import fire,os

def test(*txt):
    for i in txt:
        m=np.loadtxt(i,"U30")
        a=m[:,-1].astype("float")
        print(os.path.basename(i),*map(lambda x:np.percentile(a,x),(95,98,99)),sep="\t")

if __name__=="__main__":
    fire.Fire(test)


