## _Rana pipiens_ and the effect of pesticides
1. Transcriptome Analysis Pipeline

## Workflow:

## FastQC

1. Run FastQC

```
fastqc -o fastqc_results case_1.fastq.gz case_2.fastq.gz case_3.fastq.gz case_4.fastq.gz case_5.fastq.gz case_6.fastq.gz control_1.fastq.gz control_2.fastq.gz control_3.fastq.gz control_4.fastq.gz control_5.fastq.gz control_6.fastq.gz
```

2. Run MultiQC to Aggregate the results

```
cd fastqc_results
multiqc .
```
> Look at the MultiQC report and see if the adapters are trimmed and do a overall QC inspection.

## Index the transcriptome

## Expression Quantification

This shell script will loop through each sample and will quantify the expression
```
#!/bin/bash
for fn in ../Rpip_C{1..6};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i reference_transcriptome_index -l A ${fn}/${samp}_trim.fastq -p 8 -o quants/${samp}_quant
done

for fn in ../Rpip_T{1..6};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i reference_transcriptome_index -l A ${fn}/${samp}_trim.fastq -p 8 -o quants/${samp}_quant
done 

salmon quantmerge --quants Rpip_C1_quant Rpip_C4_quant Rpip_T1_quant Rpip_T4_quant Rpip_C2_quant Rpip_C5_quant Rpip_T2_quant Rpip_T5_quant Rpip_C3_quant Rpip_C6_quant Rpip_T3_quant Rpip_T6_quant --output merged_expression
```

## Differential Expression Analysis:

This R script will look at the correlation between samples of a treatment and plot the differential expression.

```
require(reshape2)
require(DESeq2)
expression_f <- read.table("merged_expression.tsv", header = TRUE, row.names=1)
expDesign <- data.frame(row.names = colnames(expression_f),condition = c("control","control","test","test","control","control","test","test","control","control","test","test"))

```

