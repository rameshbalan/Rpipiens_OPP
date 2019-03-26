#!/usr/bin/python
from Bio import SeqIO

de_transcripts = []

with open("DE_transcripts.txt","r") as det_file:

	for line in det_file:

		de_transcripts.append(line.rstrip())

for sequence in SeqIO.parse("reference_transcriptome.fasta", "fasta"):

	if sequence.id in de_transcripts:

		print ">"+sequence.id
		print sequence.seq
