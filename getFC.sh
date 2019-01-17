#!/bin/bash

## script to generate fold-enrichment tracks from macs2 output

prefix=$1

Sample=$(basename "$prefix")
Sample=${Sample%%.*}

sizes=$2

macs2 bdgcmp -t ${prefix}_treat_pileup.bdg -c ${prefix}_control_lambda.bdg \
      --o-prefix ${prefix} -m FE

slopBed -i ${prefix}_FE.bdg -g $sizes -b 0 | bedClip stdin $sizes ${prefix}_fc.bedgraph
rm -f ${prefix}_FE.bdg

sort -k1,1 -k2,2n ${prefix}_fc.bedgraph > ${prefix}_fc.srt.bedgraph

bedGraphToBigWig  ${prefix}_fc.srt.bedgraph $sizes ${prefix}_fc.bw

# mv  ${prefix}_fc.srt.bedgraph ${prefix}_fc.srt.bg
#bzip ${prefix}_fc.srt.bg

if [  -f ${prefix}_fc.bw ]
then
    rm ${prefix}_control_lambda.bdg
    rm ${prefix}_treat_pileup.bdg
    rm ${prefix}_fc.bedgraph
fi


