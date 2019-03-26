#!/usr/bin/bash

makeblastdb -in Xenopus_laevis_genome.fna -dbtype nucl
makeblastdb -in Xenopus_tropicalis_genome.fna -dbtype nucl

blastn -query DE_transcripts_sequence.fasta -db Xenopus_laevis_genome.fna -out DE_transcripts_vs_X_laevis.tsv -outfmt "6 qseqid sseqid evalue" -max_target_seqs 1
blastn -query DE_transcripts_sequence.fasta -db Xenopus_tropicalis_genome.fna -out DE_transcripts_vs_X_tropicalis.tsv -outfmt "6 qseqid sseqid evalue" -max_target_seqs 1