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

```
salmon index -t reference_transcriptome.fasta -i reference_transcriptome_index
```

## Expression Quantification

This shell script will loop through each sample and will quantify the expression
```
#!/bin/bash
for fn in ../Rpip_C{1..6};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i reference_transcriptome_index -l A --gcBias -r ${fn}/${samp}_trim.fastq -p 8 -o quants/${samp}_quant
done

for fn in ../Rpip_T{1..6};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i reference_transcriptome_index -l A --gcBias -r ${fn}/${samp}_trim.fastq -p 8 -o quants/${samp}_quant
done 

salmon quantmerge --quants Rpip_C1_quant Rpip_C4_quant Rpip_T1_quant Rpip_T4_quant Rpip_C2_quant Rpip_C5_quant Rpip_T2_quant Rpip_T5_quant Rpip_C3_quant Rpip_C6_quant Rpip_T3_quant Rpip_T6_quant --output merged_expression
```

## Differential Expression Analysis:
//***** In progress
*****//
This R script will look at the correlation between samples of a treatment and plot the differential expression.

```
library(DESeq2)
library(tximport)
library(rjson)
library(readr)
samples <- read.table("samples.txt", header=TRUE)
files <- file.path("/home/hdd/4/rna_rpipiens/quants", samples$sample, "quant.sf")
# txOut argument avoid genelevelsummary since gene level information is not available.
txi <- tximport(files, type="salmon",txOut=TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
ddsTxi
dds <- DESeq(ddsTxi)
res <- results(dds)
res
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_test_vs_control", type="apeglm")
resLFC
summary(res)
#How many adjusted p-values were less than 0.0.05?
sum(res$padj < 0.05, na.rm=TRUE)
```

