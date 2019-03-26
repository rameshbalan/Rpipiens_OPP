#!/usr/bin/python3


xt_ortholog_list = []
xt_ortholog_dict = {}
transcript_ann = {}


#Reading the file with has significant differential expression:
with open("DE_transcripts_vs_X_tropicalis.tsv","r") as xt_ortholog_file:

	for line in xt_ortholog_file:

		split_line = line.rstrip().split("\t")

		xt_ortholog_list.append(split_line[1])

		xt_ortholog_dict[split_line[1]] = split_line[0]


#Adding annotation based on the gff file
with open("Xenopus_tropicalis.gff","r") as gff_file:

	for line in gff_file:

		if not line.startswith("#"):

			split_line = line.split("\t")

			if split_line[2] == "mRNA":

				split_col9 = split_line[8].split(";")

				transcript_id = split_col9[-1].rstrip().split("=")

				if transcript_id[1] in xt_ortholog_list:

					print(transcript_id[1])

					transcript_ann[transcript_id[1].rstrip()] = split_line[8].split("product=")[1].split(";")[0]

#Writing the transcript ID and annotation to a new file
with open("R_pipiens_DE_X_tropicalis_orthologs.txt","w") as outfile:

	for key, value in transcript_ann.items():

		outfile.write("%s\t%s\t%s\n" % (xt_ortholog_dict[key],key,value))