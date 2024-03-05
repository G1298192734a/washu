#!/usr/bin/env python3
#-*- coding:utf-8 -*-
import fire
import glob,os

class Washu(object):
    def __init__(self,ref,species=False,url_dir="https://raw.githubusercontent.com/G1298192734a/washu/master"):
        self._r=ref
        self._s=species if species else ref
        self._u=url_dir
        self._p=lambda x,y:open(y,"w").write(x)

    def annot_tracks(self,outfile="annotationTracks.json"):
        output="""
{
    "Ruler": [
        {
            "type": "ruler",
            "label": "Ruler",
            "name": "Ruler"
        }
    ],
    "Genes": [
        {
            "name": "gene",
            "label": "%s genes",
            "filetype": "geneAnnotation"
        }
    ]
}
        """ % self._s
        self._p(output,outfile)

    def chromsize(self,size,outfile="chromSize.json"):
        csizes=list()
#        csizes=['    {"chr":"%s","size":"%s"}' % i.split() for i in open(size)]
        for i in open(size):
            chrom,size=i.split()
            csizes.append(f'    {{"chr":"{chrom}","size":"{size}"}}')
        with open(outfile,"w") as OUT:
            print("[",",\n".join(csizes),"]",sep="\n",file=OUT)

    def genome_js(self,project,jsons=".",outfile=False):
        jsons="\n".join((
            f'import {os.path.basename(i)[:-5]} from "{i}";'
            for i in glob.glob(os.path.join(jsons,"*json"))
        ))
        twobit=os.path.join(self._u,project,"genome.2bit")
        gene_annot=os.path.join(self._u,project,"genome.refbed")
        gc=os.path.join(self._u,project,"genome.gc.bigwig")
        species=self._s.replace("_","").upper()
        if not outfile:outfile=f"{self._r}.js"
        output=f"""
import Chromosome from "../Chromosome";
import Genome from "../Genome";
import TrackModel from "../../TrackModel";
{jsons}

const allSize = chromSize.map(genom => new Chromosome(genom.chr, genom.size));
const genome = new Genome("{self._r}", allSize);

const navContext = genome.makeNavContext();
const defaultRegion = navContext.parse("chr1:1000000-2000000");
const defaultTracks = [
    new TrackModel({{
        type: "ruler",
        name: "Ruler",
    }}),
    new TrackModel({{
        type: "refbed",
        name: "gene",
        label: "{self._s} genes",
        url: "{gene_annot}",
    }}),
    new TrackModel({{
        type: "bigwig",
        name: "gc",
        label: "{self._s} gc",
        url: "{gc}",
    }}),
];

const {species} = {{
    genome: genome,
    navContext: navContext,
    cytobands:{{}},
    defaultRegion: defaultRegion,
    defaultTracks: defaultTracks,
    twoBitURL: "{twobit}",
    annotationTracks,
}};

export default {species};
        """
        self._p(output,outfile)

    def hub_json(self,project,tracks=".",outfile="hub.contig.json"):
        ext_type={
        "bw":"bigwig","bb":"bigbed","gene":"bed",
        "Heatmap.sorted":"longrange","Compartment.sorted":"bed",
        "Interaction.sorted":"longrange","TAD.sorted":"bed",
        }
        jsons=list()
        for i in os.popen(f"find {tracks} -type f -name '*gz'"):
#            prefix,ext=os.path.basename(i.rstrip()).split(".")[:2]
            name=os.path.basename(i.rstrip())[:-3]
            ext=".".join(name.rsplit(".",2)[1:])
            jsons.append(f"""
    {{
    "name": "{name}",
    "type": "{ext_type.get(ext,ext)}",
    "url": "{os.path.join(self._u,project,i.rstrip())}",
    "showOnHubLoad":"true",
    }}"""
            )
        with open(outfile,"w") as OUT:
            print("[",",\n".join(jsons),"]",sep="\n",file=OUT)

if __name__=="__main__":
    fire.Fire(Washu)


