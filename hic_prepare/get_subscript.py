#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import fire
import os

class Subscript(object):
    def __init__(self,insh,p="xhacexclu12,xhacexclu16,xhacexclu03,xhacexclu11,xhacexclu09,xhacexclu13"):
        self._insh=os.path.abspath(insh)
        self._p=p if isinstance(p,str) else ",".join(p)

    def __mkdir__(self):
        tmp=f"{self._insh}_tmp"
        if os.path.exists(tmp):return
        os.mkdir(tmp)

    def __print__(self,lines):
        tmp=f"{self._insh}_tmp"
        sbatchs=f"""#!/usr/bin/sh
#SBATCH -n 1
#SBATCH -c 1
#SBTACH --mem 10G
#SBATCH -p {self._p}
cd {os.path.dirname(tmp)}
        """
        for i in range(len(lines)):
            sub_sh=os.path.join(tmp,f"{tmp.rsplit('/',1)[-1]}{i}.sh")
            with open(sub_sh,"w") as OUT:
                print(sbatchs,*lines[i],sep="\n",file=OUT)
            os.system(f"cd {tmp} && sbatch {sub_sh}") 

    def nrows(self,n):
        self.__mkdir__()
        with open(self._insh) as SH:
            ls=[i.rstrip() for i in SH if i[0].isalpha() or i.startswith("/")]
            lines=[ls[i:i+n] for i in range(0,len(ls),n)]
            self.__print__(lines)
            
    def numsign(self):
        self.__mkdir__()
        with open(self._insh) as SH:
            ls=SH.read().splitlines()
            signs=[i for i in range(len(ls)) if ls[i].startswith("#")]
            lines=[ls[signs[i]:signs[i+1]] for i in range(len(signs)-1)]
            lines.append(ls[signs[-1]:])
            self.__print__(lines)

if __name__ == "__main__":
    fire.Fire(Subscript)


