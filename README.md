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

## _De novo_ Assembly:

1. Create a Samples file. Example:

```
case	case_rep1	case_1.fastq.gz
case	case_rep2	case_2.fastq.gz
case	case_rep3	case_3.fastq.gz
case	case_rep4	case_4.fastq.gz
case	case_rep5	case_5.fastq.gz
case	case_rep6	case_6.fastq.gz
control	control_rep1	control_1.fastq.gz
control	control_rep2	control_2.fastq.gz
control	control_rep3	control_3.fastq.gz
control	control_rep4	control_4.fastq.gz
control	control_rep5	control_5.fastq.gz
control	control_rep6	control_6.fastq.gz
```

2. Run Trinity for  assembly

```
Trinity --seqType fq --samples_file samples.txt --max_memory 128G --CPU 48
```

## Check the Quality of the _de novo_ assembly vs Published transcriptome

This step will provide statistics that will help us determine which would be our reference transcriptome between the two transcriptomes.

```
transrate Trinity.fasta
```
or  

```
$TRINITY_HOME/util/TrinityStats.pl  Trinity.fasta
```

and similarly for the published transcriptome.

```
transrate published_transcriptome.fasta
```
or  

```
$TRINITY_HOME/util/TrinityStats.pl  published_transcriptome.fasta
```
## Align reads to the reference transcriptome:

1. Map each read to the reference transcriptome and estimate the read counts.

```
$TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts reference.fasta --seqType fq --samples_file samples.txt --est_method salmon --aln_method bowtie --trinity_mode --prep_reference --output_dir salmon_outdir
```

2. Combine all of them together

```
$TRINITY_HOMEutil/abundance_estimates_to_matrix.pl --est_method salmon case1_salmon.results ... control1_salmon.results ...
```

## Differential Expression Analysis:

```
$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix combined_counts.matrix --method DESeq2 --samples_file samples.txt
```

