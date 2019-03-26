library(DESeq2)
library(tximport)
library(rjson)
library(readr)
library(ashr)
library(vsn)
source("change_PCs.R")

samples <- read.table("samples.txt", header=TRUE)
files <- file.path("/home/hdd/4/rna_rpipiens/quants", samples$sample, "quant.sf")
# txOut argument avoid genelevel summary since gene level information is not available.
txi <- tximport(files, type="salmon",txOut=TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition + dev_stage)

## Differential Expression Function
dds <- DESeq(ddsTxi)
sum(res$padj < 0.05, na.rm=TRUE)
## Filter Reads with less than 10 reads across 12 samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Differential Expression Function
dds <- DESeq(dds)

# getting the results from the DE seq function
res <- results(dds, contrast = c("condition","test","control"))
resOrdered <- res[order(res$pvalue),]


#write the results to a csv
write.csv(as.data.frame(resOrdered), 
          file="condition_results.csv")

# Data Transformation
vsd <- vst(dds, blind=FALSE)

# getting the results from the DE seq function
head(assay(vsd))

# pca plot
plotPCA(vsd, intgroup="condition")

vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
mat <- limma::removeBatchEffect(mat,vsd$dev_stage)
assay(vsd) <- mat
plotPCA(vsd, intgroup="condition")
plotPCA.balan(vsd, intgroup="condition")

meanSdPlot(assay(vsd))
summary(res)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_test_vs_control", type="apeglm")
resLFC
summary(resLFC)
#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.05, na.rm=TRUE)
plotMA(resLFC, cex = 0.8)
abline(h=c(-4,4),col="dodgerblue",lwd=2)
abline(h=c(-8,8),col="skyblue",lwd=2)
abline(h=c(-10,10),col="yellow",lwd=2)
