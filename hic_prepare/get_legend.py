import fire,os,matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

def legend(pml,prefix=False,scale=10):
    # 画三维模型图例
    # 读pml，拆开chr和rgb颜色,chr一维，rgb二维
    chr_color=[i.split()[1:] for i in open(pml) if i.startswith("set_color")]
    chrs=[i[0][1:-1].capitalize() for i in chr_color] #需要pml里命名规范cchr1,去两头
    rgbs=[i[1][1:-1].split(",") for i in chr_color]
    # rgb需要转化成小数
    rgbs=[list(map(lambda x:int(x)/255.0,i)) for i in rgbs]
    # 画图例
    width,length=[i/scale for i in (max(map(len,chrs)),len(chrs))]
    plt.figure(dpi=300,figsize=(width,length))
    plt.legend(
        handles=[
            mpatches.Patch(color=i,label=j) for i,j in zip(rgbs,chrs)
        ],loc='center'
    )
    plt.axis("off")
    # 保存
    if not prefix:
        prefix=os.path.join(os.path.dirname(pml),f"legend")
    plt.savefig(f"{prefix}.png",format="png",bbox_inches='tight')#,pad_inches=0)
    plt.savefig(f"{prefix}.pdf",format="pdf",bbox_inches='tight')#,pad_inches=0)
if __name__=="__main__":
    fire.Fire(legend)


