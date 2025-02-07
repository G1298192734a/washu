#!/usr/bin/bash

module load bedtools/2.30.0
#draw_region=500000
draw_region=1000000
cont=$2
treat=$3

if [ ! -d "PGT_draw" ]; then mkdir PGT_draw; fi
cd PGT_draw
while read -r -a line; do 
  if [ ! -d ${line[4]}_${line[5]} ]; then mkdir ${line[4]}_${line[5]}; fi
  cd ${line[4]}_${line[5]}
  if [ ! -d 00.prepare_data ]; then ln -s ../../../00.prepare_data; fi
  chr_name=${line[0]}
  matrix1=00.prepare_data/01.cool_bychr/${treat}_40k.balance.cool
  tad1=00.prepare_data/04.tad/${treat}/corrected_${chr_name}.TAD.bed
  matrix2=00.prepare_data/01.cool_bychr/${cont}_40k.balance.cool
  tad2=00.prepare_data/04.tad/${cont}/corrected_${chr_name}.TAD.bed
  minus_matrix=00.prepare_data/02.minus_matrix_cool/${cont}_vs_${treat}/${chr_name}_zscore_minus.cool
  cpt1=00.prepare_data/03.compartment/${treat}_100k_cpt.bed
  cpt2=00.prepare_data/03.compartment/${cont}_100k_cpt.bed
  insulation1=00.prepare_data/05.tad_scores/${treat}.40000.score.bw
  insulation2=00.prepare_data/05.tad_scores/${cont}.40000.score.bw
  vline=00.prepare_data/04.tad/${cont}_vs_${treat}/corrected_${chr_name}.merged.TAD.bed
  #RNA_1=00.prepare_data/06.RNA_bw/${treat}-1.bw
  if [ ! -f ${line[4]}_gene.bed ]; then printf "%s\t%s\t%s\t%s\n" "${line[0]}" "${line[1]}" "${line[2]}" "${line[4]}" > ${line[4]}_gene.bed; fi
  if [ ! -f ${line[4]}_slop.bed ]; then perl -anle '$F[1]=int(($F[1]+$F[2])/2);$F[2]=$F[1]+1;print join "\t",@F' ${line[4]}_gene.bed|slopBed -b $draw_region -i - -g ../../genome.fa.fai|perl -anle "if(\$F[1]==0){\$F[2]=\$F[1]+2*$draw_region}else{\$F[1]=\$F[2]-2*$draw_region} print join qq(\t),@F" > ${line[4]}_slop.bed; fi
#  if [ ! -f ${treat}_${chr_name}_target.arcs ]; then pairToBed -a 00.prepare_data/06.loop/${treat}_cis.bedpe -b ${line[4]}_gene.bed -type either|cut -f1-7|sort|uniq > ${treat}_${chr_name}_target.arcs; fi
  if [ ! -f ${treat}_${chr_name}_regional.arcs ]; then pairToBed -a 00.prepare_data/06.loop/${treat}_cis.bedpe -b ${line[4]}_slop.bed -type both|cut -f1-7|sort|uniq > ${treat}_${chr_name}_regional.arcs; fi
#  if [ ! -f ${cont}_${chr_name}_target.arcs ]; then pairToBed -a 00.prepare_data/06.loop/${cont}_cis.bedpe -b ${line[4]}_gene.bed -type either|cut -f1-7|sort|uniq > ${cont}_${chr_name}_target.arcs; fi
  if [ ! -f ${cont}_${chr_name}_regional.arcs ]; then pairToBed -a 00.prepare_data/06.loop/${cont}_cis.bedpe -b ${line[4]}_slop.bed -type both|cut -f1-7|sort|uniq > ${cont}_${chr_name}_regional.arcs; fi
  backslash="\\"
  tracks=$(cat <<- END_HEREDOC
[hicmatrix]
file = ${matrix1}
#depth = 200000
#depth = 500000
#depth = 400000
depth = 700000
#depth = 900000
title = ${treat}
#transform = log1p
transform = log
file_type = hic_matrix
show_masked_bins = false
#show_masked_bins = true
colormap = RdYlBu_r
height = 8

[tads]
file = ${tad1}
display = triangles
border_color = black
#border_color = red
color = none
line_width = 1.5
overlay_previous = share-y
height = 8

[spacer]
height = 0.1

[hicmatrix]
file = ${matrix2}
#depth = 200000
#depth = 500000
#depth = 400000
depth = 700000
#depth = 800000
#depth = 900000
title = ${cont}
#transform = log1p
transform = log
file_type = hic_matrix
show_masked_bins = false
#show_masked_bins = true
colormap = RdYlBu_r
height = 8
orientation = inverted

[tads]
file = ${tad2}
display = triangles
border_color = black
#border_color = red
color = none
line_width = 1.5
overlay_previous = share-y
#height = 5
height = 8

[spacer]

[hicmatrix]
file = ${minus_matrix}
#depth = 200000
#depth = 600000
#depth = 500000
depth = 800000
title = ${cont}_vs_${treat} matrix
#min_value = -1
#max_value = 1
min_value = -2
max_value = 2
transform = no
file_type = hic_matrix
colormap = RdBu_r
show_masked_bins = false
height = 10

[spacer]

[gff genes]
file = ../../genome.gtf
height = 15 
fontsize = 10
#merge_transcripts = true
prefered_name = gene_name
file_type = gtf
max_labels = 2000
style = flybase

[bed genes]
file = ${line[4]}_gene.bed
height = 1
fontsize = 15
style = flybase
color = green
title = target region

[spacer]
height = 0.5

[cpt1]
file = ${cpt1}
height = 1.5
title = ${treat} comparment A/B
min_value = -1
max_value = 1
color = red
negative_color = blue
file_type = bedgraph

[spacer]
height = 1.2

[cpt2]
file = ${cpt2}
height = 1.5
title = ${cont} comparment A/B
min_value = -1
max_value = 1
color = red
negative_color = blue
file_type = bedgraph

[spacer]
height = 0.5

[arcs1]
file = ${treat}_${chr_name}_regional.arcs
title = ${treat} contacts
height = 5
color = RdYlBu
line_width = 1.5
use_middle = true

[hlines]
height = 0.3
show_data_range = false
y_values = -0.3, 0.3
file_type = hlines

[arcs2]
file = ${cont}_${chr_name}_regional.arcs
title = ${cont} contacts
height = 5
color = RdYlGn
line_width = 1.5
use_middle = true
orientation = inverted

[spacer]

[insulation 1]
title = insulation:
  ${treat} red
  ${cont} blue
grid = true
color = red
type = line
file_type = bigwig
file = ${insulation1}
number_of_bins = 70
height = 5
nans_to_zeros = true
show_data_range = true
summary_method = max

[insulation 2]
grid = true
color = #11c5fd
type = line
file_type = bigwig
file = ${insulation2}
number_of_bins = 70
height = 5
overlay_previous = share-y
nans_to_zeros = true
show_data_range = true
summary_method = max

[vlines]
file = ${vline} 
type = vlines
[x-axis]
END_HEREDOC
)
  printf "%s" "${tracks}" > ${cont}_vs_${treat}_track_${draw_region}.ini
  python3 /work/frasergen/3D/work/shaojie/script/HiC/integrate/modify_ini.py ${cont}_vs_${treat}_track_${draw_region}.ini ${line[4]}_slop.bed
  draw_start=$(head -1 ${line[4]}_slop.bed|cut -f2)
  draw_end=$(head -1 ${line[4]}_slop.bed|cut -f3)
  region=$[$draw_end-$draw_start]
  if [ $draw_start = 0 ]; then
    draw_end=$[2*$draw_region]
  else
    draw_start=$[$draw_end-2*$draw_region]
  fi
  cmd=$(cat <<- END_DOC
#!/usr/bin/bash
/public/frasergen/3D/work/wanghao/mamba/PGT-3.9/bin/pyGenomeTracks --tracks ${cont}_vs_${treat}_track_${draw_region}.modi.ini --trackLabelFraction 0.08 --region ${chr_name}:${draw_start}-${draw_end} -o ${cont}_vs_${treat}_${line[4]}_${line[5]}.png --height 60 --width 50 --dpi 300
/public/frasergen/3D/work/wanghao/mamba/PGT-3.9/bin/pyGenomeTracks --tracks ${cont}_vs_${treat}_track_${draw_region}.modi.ini --trackLabelFraction 0.08 --region ${chr_name}:${draw_start}-${draw_end} -o ${cont}_vs_${treat}_${line[4]}_${line[5]}.pdf --height 60  --width 50 --dpi 300
#/public/frasergen/3D/work/wanghao/myconda/stripenn/bin/python /public/frasergen/backup/3d/project/Interactome/211600_mianyang/09.personal/Hi-C_integrate_plots_20230906/02.heatmap_matrix_plot/matrix_draw_subplots.py -c1 ${treat}.cool -c2 ${cont}.cool -r 20000 -c chr1 -s ${draw_start} -e ${draw_end} -s1 $treat -s2 $cont -b False
END_DOC
)
  if [ ! -f ${cont}_vs_${treat}_PGT_draw.sh ]; then printf "%s" "${cmd}" > ${cont}_vs_${treat}_PGT_draw.sh;fi
  sh ${cont}_vs_${treat}_PGT_draw.sh && echo ${cont}_vs_${treat}_${line[4]}_${line[5]} "PGT draw done"
  
  cd ..
done < "$1" 

