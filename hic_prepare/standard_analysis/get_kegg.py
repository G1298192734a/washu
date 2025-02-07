import fire,yaml
from pathlib import Path

config=yaml.load(open("/work/frasergen/3D/pipeline/Interactome/hic3_workflow/database.yaml", "r"), Loader=yaml.FullLoader)
def kegg(prefix,hic="./",species="a"):
    a={"a":"animal","p":"plant","f":"fungi"}
    s=a[species]
    hic=Path(hic).resolve()
    ref=hic.parent.joinpath("00.ref")
    hic=hic.joinpath("01.ref")
    KEGG,NR=hic.joinpath("KEGG"),hic.joinpath("NR")
    for i in (KEGG,NR):
        i.mkdir(exist_ok=True,parents=True)
    script=f"""#!/usr/bin/sh
/public/frasergen/PUB/software/python/Python-3.7.7/bin/python3 /public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/chrom_trans.py --inref {ref.joinpath(prefix)}.fasta --gene-bed {ref.joinpath(prefix)}_gene.bed --chrom-map {ref.joinpath("mapfile_hic.txt")} --outfa {hic.joinpath("genome.fa")} --outgene {hic.joinpath("genome_gene.bed")}
/public/frasergen/PUB/software/samtools/samtools-1.12/bin/samtools faidx {hic.joinpath("genome.fa")}
/public/frasergen/PUB/software/bwa/bwa-0.7.17/bwa index {hic.joinpath("genome.fa")}
/public/frasergen/PUB/software/UCSCTOOLS/faSize -detailed -tab {hic.joinpath("genome.fa")} > {hic.joinpath("genome.chrsize")}
    """
    open(hic.joinpath("work.sh"),"w").write(script)
    script=f"""#!/usr/bin/sh
module load diamond/2.0.7 ;
diamond blastx --db {config[s+"_kegg"]} --query {ref.joinpath(prefix)}.mrna.fasta -c1 --threads 8 --outfmt 6 --evalue 1e-5 --more-sensitive --max-target-seqs 15 --out {KEGG.joinpath("genome.gene.ko")}
mypath=$PATH ;
export PATH=/public/frasergen/RNA/pipeline/software/Anaconda2/anaconda2/bin:$PATH ;
python /public/frasergen/RNA/pipeline/smallRNA_Pipeline_V1.0/lib/prepare/lib/kobas-3.0/scripts/annotate.py -i {KEGG.joinpath("genome.gene.ko")} -t blastout:tab -s ko -k /public/frasergen/RNA/pipeline/smallRNA_Pipeline_V1.0/lib/prepare/lib/kobas-3.0 -q /public/frasergen/RNA/pipeline/smallRNA_Pipeline_V1.0/lib/prepare/lib/sqlite3/ -o {KEGG.joinpath("genome.gene.annot")}
export PATH=/public/frasergen/PUB/software/python/Python-2.7.9/bin:$mypath ;
/public/frasergen/PUB/software/python/Python-3.7.7/bin/python3 /public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/kegg2go.py --kokeg /public/frasergen/PUB/database/KEGG/20191023/ko00001.keg --gene2mark {ref.joinpath(prefix)}.gene2mark --kegg-annot {KEGG.joinpath("genome.gene.annot")} --outpfix {KEGG.joinpath("genome")}
    """
    open(KEGG.joinpath("work.sh"),"w").write(script)
    script=f"""#!/usr/bin/sh
module load diamond/2.0.7 ;
diamond blastx --db {config[s+"_nr"]} --query {ref.joinpath(prefix)}.mrna.fasta --threads 8 -c1 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle --more-sensitive --evalue 1e-5 --max-target-seqs 30 --salltitles --out {NR.joinpath("genome.gene.nraln")}
/public/frasergen/PUB/software/python/Python-3.7.7/bin/python3 /public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/nr2go.py --nr-align {NR.joinpath("genome.gene.nraln")} --access2go /public/frasergen/PUB/database/GO/2021-09-01/accession2go --gotab /public/frasergen/PUB/database/GO/2021-09-01/go.tab --gene2mark {ref.joinpath(prefix)}.gene2mark --outpfix {NR.joinpath("genome")}
/public/frasergen/PUB/software/perl/perl-5.34.0/bin/perl /public/frasergen/3D/work/shaojie/pipeline/ThreeDgenome/03.Interact/ICE/databases/enrich/bin/go_classification_forall.pl {NR.joinpath("genome.wego")} {NR.joinpath("genome")}
/public/frasergen/PUB/software/perl/perl-5.34.0/bin/perl /public/frasergen/3D/work/shaojie/pipeline/ThreeDgenome/03.Interact/ICE/databases/enrich/bin/get_GO_anno.pl -go {NR.joinpath("genome.GO_classification_all.xls")} -o {NR.joinpath("genome.ACC2GO")}
    """
    open(NR.joinpath("work.sh"),"w").write(script)

if __name__=="__main__":
    fire.Fire(kegg)

