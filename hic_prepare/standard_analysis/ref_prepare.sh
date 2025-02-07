#!/bin/bash
if [ $# -eq 0 ]; then 
	echo "PROGRAM| reference prepare for Hi-C3.0 pipeline "
    echo "  USAGE| ref_prepare.sh <raw_ref> <raw_gff> <chrom_map> <meme_db> <outpfix>"
	echo "         <meme_db> used to scan motif in whole genome by fimo, "
	echo "         eg. /public/frasergen/PUB/software/meme/meme-5.5.5/database/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
	echo 
	echo "lurui@frasergen.com"
	exit 1
fi

raw_fasta=$1
raw_gff=$2
chrom_map=$3
outpfix=$4


python3 /public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/hic3_util.py fa_clean --infa $raw_fasta --chrom-map $chrom_map --outpfix $outpfix
echo "fa clean: ok !"
python3 /public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/hic3_util.py gff_clean --ingff $raw_gff --chrom-map $chrom_map --outpfix $outpfix
echo "gff clean: ok !"
python3 /public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/hic3_util.py gff_gene --ingff ${outpfix}.gff --outpfix $outpfix
echo "gff convert to genome bed: ok !"
gffread -T ${outpfix}.gff -o ${outpfix}.gtf
echo "gff convert to gtf: ok !"
python3 /public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/hic3_util.py gtf_gene2mark --ingtf ${outpfix}.gtf --outpfix ${outpfix}
echo "gene2mark generate: ok !"
gffread ${outpfix}.gtf -g ${outpfix}.fasta -x ${outpfix}.mrna.fasta
echo "get transcript fasta file: ok !"
#python3 /public/frasergen/3D/pipeline/Interactome/hic3_workflow/scripts/parallel_fimo.py --infasta ${outpfix}.fasta --meme_db $meme_db --outdir $(dirname $outpfix)
#echo "fimo scan whole genome: ok!"
