'''This program will convert a maf transcriptome alignment to fasta format.
It will filter out sequences with less than 60 nucleotides (shortest known 
protein is 20 residues long).

Note: it will not work for genome alingments since it will include all 
sequences in the output. '''

import argparse
from mafTrans import convert

def main():
	parser = argparse.ArgumentParser(description = "This program will convert \
a maf transcriptome alignment to fasta format. It will filter out sequences \
with less than 60 nucleotides. Note: it will not work for genome \
alingments since it will include all sequences in the output.")
	parser.add_argument("i", help = "Path to input maf file.")
	args = parser.parse_args()
	outfile = args.i.replace(".maf", ".fasta")
	convert(args.i, outfile)

if __name__ == "__main__":
	main()
