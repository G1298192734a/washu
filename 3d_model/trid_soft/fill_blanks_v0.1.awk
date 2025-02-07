#!/usr/bin/awk
#v0.1 add chrs parse!!! input res and fasta.fai to parse per chr lines!!!!
#usage: awk -v res=${res} -f fill_blanks.awk ${genome.fa.fai} ${in.heatmap} > ${out.heamap} 
#usage i.e. awk -v res=400000 -f fill_blanks_v0.1.awk genome.fa.fai heatmap_chr21_trimmed > heatmap_chr21_trimmed_filled
# simple fill zeros with above lines(same as columns)
function round(x){
    return x==int(x)?x:int(x)+1
}
#parse blank line in bins pos(line!!)
function parse_blank_chr(chr_lines,chrl_counts,in_line){
    for(z = 1;z<=chrl_counts;z++){
        if(in_line <= chr_lines[z] && in_line > chr_lines[z-1]){
            return z;
            break;
        }
    }
}
BEGIN{
    #add chroms[0] for blank calculation!!!
    chroms[0]=0
}
NR==FNR{

    #now start to parse chr bin pos!!
    chr_count++
    chroms[chr_count]= round($2/res) + line_cash
    line_cash = chroms[chr_count]
    next;
}
#END{for(i in chroms){
#   print i , chroms[i]}
#   print chr_count
#   #test blank test
#   print parse_blank_chr(chroms,chr_count, 200)
#   print parse_blank_chr(chroms,chr_count, 1)
#   print parse_blank_chr(chroms,chr_count, 488)
#   print parse_blank_chr(chroms,chr_count, 489)
#   print parse_blank_chr(chroms,chr_count, 4726)
#   print parse_blank_chr(chroms,chr_count, 4727)
#   print parse_blank_chr(chroms,chr_count, 4728)
#   print parse_blank_chr(chroms,chr_count, 6819)
#}

#parse raw txt matrix
{
    for(i = 1;i <= NF; i++){
        $i == 0 ? zeros++ : "";
        #store matrix value
        matrix[FNR][i] = $i;
    }
    #store each chrs blanks as 2D array
    if(zeros/NF > 0.90){
        blanks[parse_blank_chr(chroms,chr_count,FNR)][FNR];
        #print("first test\t" FNR)
    }
    zeros = 0;
    matrix_scale = NF;
}
END{
    #determin blank lines change to which above lines ### edited to parse per chr info
    for(chrs=1;chrs<=chr_count;chrs++){
    #test whether this chr have blank lines!!!
        if(isarray(blanks[chrs])){
            #print("####\t" chrs)
            for(a in blanks[chrs]){
                #recursive search for not blank lines!! caution for i(line) line numbers!!! 1st search above, if first line blank or in before chr then search down 
                for(i = a; i in blanks[chrs]; i--){
                     #test whether hit line 0 or last chr...
                     if(i - 1 != 0 && i - 1 != chroms[chrs - 1] ){
                         #array store blank lines changed to not empty lines info
                         blanks[chrs][a] = i - 1
                     #if hit last chr!!!! then search donw use i++ !!
                     }else{
                         for(i = a; i in blanks[chrs]; i++){
                             blanks[chrs][a] = i + 1
                         }
                         #i++ to exit recursive loops
                         i++
                     }
                }
            }
        }
    }
    #####!!!!now test blank lines changed to info !!!!####
    #for(chrs=1;chrs<=chr_count;chrs++){
    #    if(isarray(blanks[chrs])){
    #        for(a in blanks[chrs]){
    #            print chrs "\t" a "\t" blanks[chrs][a]
    #        }
    #    }
    #}

#change values in blank lines per chrs
for(chrs=1;chrs<=chr_count;chrs++){ 
    #test blanks array not empty
    if(isarray(blanks[chrs])){
        for(a in blanks[chrs]){
            #change rows and columns values
            for(i = 1; i <= matrix_scale; i++){
                #change rows
                matrix[a][i] = matrix[blanks[chrs][a]][i];
                #change columns
                matrix[i][a] = matrix[i][blanks[chrs][a]];
            }
        
        }
    }
}
#print out matrix
for(i = 1; i <= matrix_scale; i++){
    for(j = 1; j < matrix_scale; j++){
        printf("%s ", matrix[i][j]);
    }
    #print last columns wiht "\n"
    printf("%s\n", matrix[i][matrix_scale])
}

}
