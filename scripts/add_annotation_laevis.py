#!/usr/bin/python3

#Declaring list and dictionary
xl_ortholog_list = []
xl_ortholog_dict = {}
xt_ortholog_list = []
xt_ortholog_dict = {}
transcript_ann = {}

#Reading the file with has significant differential expression:
with open("DE_transcripts_vs_X_laevis.tsv","r") as xl_ortholog_file:

	for line in xl_ortholog_file:

		split_line = line.rstrip().split("\t")

		xl_ortholog_list.append(split_line[1])

		xl_ortholog_dict[split_line[1]] = split_line[0]

#Adding annotation based on the gff file
with open("Xenopus_laevis.gff","r") as gff_file:

	for line in gff_file:

		if not line.startswith("#"):

			split_line = line.split("\t")

			if split_line[2] == "mRNA":

				split_col9 = split_line[8].split(";")

				transcript_id = split_col9[-1].rstrip().split("=")

				#print(transcript_id)

				if transcript_id[1] in xl_ortholog_list:

					transcript_ann[transcript_id[1].rstrip()] = split_line[8].split("product=")[1].split(";")[0]

#Writing the transcript ID and annotation to a new file
with open("R_pipiens_DE_X_laevis_orthologs.txt","w") as outfile:

	for key, value in transcript_ann.items():

		outfile.write("%s\t%s\t%s\n" % (xl_ortholog_dict[key],key,value))
