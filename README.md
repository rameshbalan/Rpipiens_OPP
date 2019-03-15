# _Rana pipiens_ Transcriptome Analysis Pipeline:

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
#Importing All the libraries
library(DESeq2)
library(tximport)
library(rjson)
library(readr)
library(ashr)

#Reading in the samples file which describes each sample as control/test.
samples <- read.table("samples.txt", header=TRUE)

#Reading in all the quantification files from all the 12 samples
files <- file.path("/home/hdd/4/rna_rpipiens/quants", samples$sample, "quant.sf")
# txOut argument avoid genelevelsummary since gene level information is not available.
txi <- tximport(files, type="salmon",txOut=TRUE)
#creating a DE dataset object using the samples file and all the quantification files
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
ddsTxi

# Running Differential Expression
dds <- DESeq(ddsTxi)
# Filtering genes with less than or equal to 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Running Differential Expression (Again...)
dds <- DESeq(dds)

# Getting the results from the dataset object 
res <- results(dds)
res

#Ordering the file based on the pvalue
resOrdered <- res[order(res$pvalue),]

# Provides a summary on the number of genes upregulated or downregulated
summary(res)
resultsNames(dds)

#Shrinks/Reduces the log fold change(LFC) based on the the control vs test group using apeglm model.  
resLFC <- lfcShrink(dds, coef="condition_test_vs_control", type="apeglm")
resLFC
summary(resLFC)

#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

# Getting the results from the dataset object. The default cutoff is 0.1 so changing alpha to 0.05. 
res05 <- results(dds, alpha=0.05, name = "condition_test_vs_control")
summary(res05)

# Get the number of genes differentially expressed
sum(res05$padj < 0.05, na.rm=TRUE)

# Plot the fold change 
plotMA(resLFC, ylim=c(-10,10))
#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]

```

## Preliminary Results

1. Number of significantly Upregulated and Downregulated Genes.

```
out of 142585 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)     : 358, 0.25% 
LFC < 0 (down)   : 351, 0.25% 
outliers [1]     : 22711, 16% 
low counts [2]   : 66668, 47% 
(mean count < 5)
```

2. Total Number of Differentially Expressed Genes : `709`
