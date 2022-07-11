#!/bin/env python3


import os, sys, argparse


def parseFasta(filename):
	# Open file for reading
	header = None
	seq = ""
	for line in open(filename):
		# Check if new entry
		if line.startswith('>'):
			# New entry found, yield old, reset
			if header is not None:
				yield (header, seq)
			header = line
			seq = ""
		else:
			# Sequence line
			seq += line

	# Yield last line
	if header is not None:
		yield (header, seq)


def main():
	parser = argparse.ArgumentParser(description="Rename assembly and gff file based a nomenclature file.", allow_abbrev=False)
	parser.add_argument("--nomenclature_file", type=str, required=True, 
		help="A tab deliminated file containing 2 columns - 'old-name\tnew-name'")
	parser.add_argument("--fasta", type=str, default="",
		help="The input fasta file (default: no)")
	parser.add_argument("--output_fasta", type=str, default="",
		help="The output fasta file (default: no)")
	parser.add_argument("--gff_file", type = str, default = "",
		help = "Annotation file in a gff format. Default: No")
	parser.add_argument("--output_gff_file", type = str, default = "",
		help = "Updated annotation file in a gff format. Default: No")

	args = parser.parse_args()
	# Open outputfiles
	outfile_fa = open(args.output_fasta,"w")
	outfile_gff = open(args.output_gff_file,"w")

	# args.nomenclature_file is a tab deliminated file with 2 columnes 'old-name\tnew-name' 
	# The <old-name> in input fasta will be replaced by <new-name>

	# Read and store the information from the nomenclature file
	dictionary_nomenclature = {}
	for line in open(args.nomenclature_file):
		old_name, new_name = line.rstrip('\n').split("\t")
		dictionary_nomenclature[old_name] = new_name

	# Read the input fasta file and create a new fasta file containing the new names
#	for line in open(args.fasta):
#		if line.startswith(">"):
#			name_fasta = line[1:].strip()
#			print("fasta header:" + name_fasta)
#			if name_fasta in dictionary_nomenclature:
#				line = '>' + dictionary_nomenclature[name_fasta] + '\n'
#		outfile_fa.write(line)

	for header, seq in parseFasta(args.fasta):
		# Get name
		name_fasta = header[1:].strip()
		if name_fasta in dictionary_nomenclature:
			line = '>' + dictionary_nomenclature[name_fasta] + '\n'
			outfile_fa.write(line)
			outfile_fa.write(seq)


	# Read the gff annotation file and update the chromosome names
	for line in open(args.gff_file):
		liftoff_chr_name, origin_annotation, features, start, end, junk1, direction, junk2, info = line.rstrip('\n').split("\t")
		if liftoff_chr_name in dictionary_nomenclature:
			liftoff_new_name = dictionary_nomenclature[liftoff_chr_name]
			outfile_gff.write(liftoff_new_name + "\t" + origin_annotation + "\t" + features + "\t" + start + "\t" + end + "\t" + junk1 + "\t" + direction + "\t" + junk2 + "\t" + info + "\n")

if __name__ == "__main__":
	main()