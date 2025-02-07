#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import fire
from datetime import datetime as dt,timedelta
from pandas import read_excel
#import pandas as pd
month=dt.now().strftime("%Y-%m")
today=dt.now()
monday=dt.strftime(today-timedelta(today.weekday()),"%Y-%m-%d")
friday=dt.strftime(today+timedelta(4-today.weekday()),"%Y-%m-%d")

class Summary(object):
    def week(self,inxls,outxls="week_summary.xlsx"):
        df=read_excel(inxls,usecols=[
            "合同编号","合同名称","开始时间","标题","产品类型",
        ])
        df=df[[monday<=i<=friday for i in df["开始时间"]]]
        df["任务"]=df["合同编号"]+" "+df["合同名称"]+" "+df["产品类型"]
        df["情况"]="100%"
        df["内容"]=df["标题"]
        df.to_excel(outxls,sheet_name="sheet1",index=False)
    
    def month(self,inxls,month=month,outxls="month_summary.xlsx"):
        df=read_excel(inxls,usecols=[
            "合同编号","合同名称","开始时间","任务单状态","标题","产品类型","实际完成时间",
        ])
        #print([i for i in df["开始时间"]]);exit()
        ###df=df[[i.startswith(month) for i in df["开始时间"]]]
        hic=df[[str(i).endswith("Hi-C分析") for i in df["产品类型"]]]
        hics=[",".join(hic.loc[i][["合同编号","合同名称","产品类型","任务单状态"]]) for i in hic.index]
        df["成果"]=df["汇总"]=""
        df["成果"].loc[1]="\r\n".join([(f"{i+1}、{hics[i]}") for i in range(len(hics))])
        df["汇总"].loc[1]=f"{len(hics)}个分析，完成"
        after=df[[not str(i).endswith("Hi-C分析") for i in df["产品类型"]]]
        afters=[",".join(after.loc[i][["合同编号","合同名称","产品类型","任务单状态"]]) for i in after.index]
        df["成果"].loc[2]="\r\n".join([(f"{i+1}、{afters[i]}") for i in range(len(afters))])
        df["汇总"].loc[2]=f"{len(afters)}次售后，完成"
        df=df.reindex(columns=[
            "合同编号","合同名称","产品类型","标题","开始时间","实际完成时间","成果","汇总",
        ])
        df["实际完成时间"]=[str(i)[:10].replace("-","/") for i in df["实际完成时间"]]
        df["开始时间"]=df["开始时间"].str.replace("-","/")
        df.insert(4,"样本数","*")
        df.insert(5,"工时","4")
        df.insert(8,"填表人","邵杰")
        df.to_excel(outxls,sheet_name="sheet1",index=False)

if __name__=="__main__":
    fire.Fire(Summary)


