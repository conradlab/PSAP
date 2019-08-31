#!/usr/bin/env python3
"""
This script take the ANNOVAR annoated files (either potentiomes or exome data) and split one fusion gene row into two, annotated each with one of the fusion gene.  

Usage:
	python3 fusion_gene_split.py <input file>

	Arguments:
	<Ainput file> is the annotated text files after ANNOVAR annotation, need to provide the path to the file
"""
# import modules
import sys
import os
#Make sure # of arguments is correct
if len(sys.argv) != 2:
	sys.exit( __doc__)
#Extract file name and directory
input = sys.argv[1]
file_path = os.path.dirname(input)
print("file_path = "+ file_path)
file_name = os.path.splitext(input)[0] #file_name = os.path.basename(input).split(".")[0]
print("file_name = "+ file_name)

## Read in annotated file
data_file = open(input,"r")
# read the header line and find the column index for Gene, this information is used later for duplicate variant and slipt gene
header = data_file.readline().rstrip("\n").split("\t") #read the current line, strip the end "\n", and split by "\t"
Gene_ind = [i for i, j in enumerate(header) if "Gene.wgEncodeGencodeBasic" in j][0] #find the index of column for Gene annotation
## Open output file
out_file = open(file_name + "_dup.txt","w")
# write out the same header as the input file
out_file.write("\t".join(header) + "\n")
## Read through lines
for lines in data_file:
	# Split lines
	line = lines.rstrip("\n").split("\t")
	# extracts gene info
	Gene = line[Gene_ind]
	# two situations: fusion gene or not
	if ";" not in Gene:	#not fusion gene
		out_file.write(lines)
	else: #fusion gene
		Genes = Gene.split(";") # a list of all the genes in the fusion gene, sometimes > 2 genes
		for i in range(len(Genes)):
			line_new = line
			line_new[Gene_ind] = Gene.split(";")[i]
			out_file.write("\t".join(line_new) + "\n")
out_file.close()
data_file.close()