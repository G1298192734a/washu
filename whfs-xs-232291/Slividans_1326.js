
import Chromosome from "../Chromosome";
import Genome from "../Genome";
import TrackModel from "../../TrackModel";
import annotationTracks from "./annotationTracks.json";
import chromSize from "./chromSize.json";

const allSize = chromSize.map(genom => new Chromosome(genom.chr, genom.size));
const genome = new Genome("Slividans_1326", allSize);

const navContext = genome.makeNavContext();
const defaultRegion = navContext.parse("chr1:1000000-2000000");
const defaultTracks = [
    new TrackModel({
        type: "ruler",
        name: "Ruler",
    }),
    new TrackModel({
        type: "refbed",
        name: "gene",
        label: "Slividans_1326 genes",
        url: "https://raw.githubusercontent.com/G1298192734a/washu/master/whfs-xs-232291/genome.refbed.gz",
    }),
    new TrackModel({
        type: "bigwig",
        name: "gc",
        label: "Slividans_1326 gc",
        url: "https://raw.githubusercontent.com/G1298192734a/washu/master/whfs-xs-232291/genome.gc.bigwig",
    }),
];

const SLIVIDANS1326 = {
    genome: genome,
    navContext: navContext,
    cytobands:{},
    defaultRegion: defaultRegion,
    defaultTracks: defaultTracks,
    twoBitURL: "https://raw.githubusercontent.com/G1298192734a/washu/master/whfs-xs-232291/genome.2bit",
    annotationTracks,
};

export default SLIVIDANS1326;
        
