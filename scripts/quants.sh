#!/bin/bash
for fn in Rpip_C{1..6};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i reference_transcriptome_index -l A --gcBias -r ${samp}_trim.fastq -p 16 -o quants/${samp}_quant
done

for fn in Rpip_T{1..6};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i reference_transcriptome_index -l A --gcBias -r ${samp}_trim.fastq -p 16 -o quants/${samp}_quant
done 
