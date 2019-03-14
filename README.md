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

## Align reads to the reference transcriptome:

1. Map each read to the reference transcriptome and estimate the read counts.

```
bwa index reference_transcriptome.fasta
```

```
#!/usr/bin/bash

bwa aln reference_transcriptome.fasta Rpip_C1_trim.fastq > Rpip_C1_trim.sai
bwa samse reference_transcriptome.fasta Rpip_C1_trim.sai Rpip_C1_trim.fastq > Rpip_C1_trim.sam

bwa aln reference_transcriptome.fasta Rpip_C2_trim.fastq > Rpip_C2_trim.sai
bwa samse reference_transcriptome.fasta Rpip_C2_trim.sai Rpip_C2_trim.fastq > Rpip_C2_trim.sam

bwa aln reference_transcriptome.fasta Rpip_C3_trim.fastq > Rpip_C3_trim.sai
bwa samse reference_transcriptome.fasta Rpip_C3_trim.sai Rpip_C3_trim.fastq > Rpip_C3_trim.sam

bwa aln reference_transcriptome.fasta Rpip_C4_trim.fastq > Rpip_C4_trim.sai
bwa samse reference_transcriptome.fasta Rpip_C4_trim.sai Rpip_C4_trim.fastq > Rpip_C4_trim.sam

bwa aln reference_transcriptome.fasta Rpip_C5_trim.fastq > Rpip_C5_trim.sai
bwa samse reference_transcriptome.fasta Rpip_C5_trim.sai Rpip_C5_trim.fastq > Rpip_C5_trim.sam

bwa aln reference_transcriptome.fasta Rpip_C6_trim.fastq > Rpip_C6_trim.sai
bwa samse reference_transcriptome.fasta Rpip_C6_trim.sai Rpip_C6_trim.fastq > Rpip_C6_trim.sam

bwa aln reference_transcriptome.fasta Rpip_T1_trim.fastq > Rpip_T1_trim.sai
bwa samse reference_transcriptome.fasta Rpip_T1_trim.sai Rpip_T1_trim.fastq > Rpip_T1_trim.sam

bwa aln reference_transcriptome.fasta Rpip_T2_trim.fastq > Rpip_T2_trim.sai
bwa samse reference_transcriptome.fasta Rpip_T2_trim.sai Rpip_T2_trim.fastq > Rpip_T2_trim.sam

bwa aln reference_transcriptome.fasta Rpip_T3_trim.fastq > Rpip_T3_trim.sai
bwa samse reference_transcriptome.fasta Rpip_T3_trim.sai Rpip_T3_trim.fastq > Rpip_T3_trim.sam

bwa aln reference_transcriptome.fasta Rpip_T4_trim.fastq > Rpip_T4_trim.sai
bwa samse reference_transcriptome.fasta Rpip_T4_trim.sai Rpip_T4_trim.fastq > Rpip_T4_trim.sam

bwa aln reference_transcriptome.fasta Rpip_T5_trim.fastq > Rpip_T5_trim.sai
bwa samse reference_transcriptome.fasta Rpip_T5_trim.sai Rpip_T5_trim.fastq > Rpip_T5_trim.sam

bwa aln reference_transcriptome.fasta Rpip_T6_trim.fastq > Rpip_T6_trim.sai
bwa samse reference_transcriptome.fasta Rpip_T6_trim.sai Rpip_T6_trim.fastq > Rpip_T6_trim.sam
```


## Differential Expression Analysis:

```
$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix combined_counts.matrix --method DESeq2 --samples_file samples.txt
```

