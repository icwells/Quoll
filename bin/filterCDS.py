'''This program will sort TransDecoder output to identify open reading frames 
in a nucleotide alignment and remove genes with insufficient coverage.'''

import argparse
from findCDS import *
			 
def main():
	parser = argparse.ArgumentParser(description = "This program will \
identify open reading frames in a nucleotide alignment and remove genes with \
insufficient coverage.")
	parser.add_argument("i", help = "Path to input file.")
	parser.add_argument("-p", default = 0.9, type = float,
help = "Minimum nucleotide content for each coding sequence (default = 0.9).") 
	args = parser.parse_args()
	outfile = args.i[:args.i.find(".")] + ".ORFs.fasta"
	index = indexAlignment(args.i)
	cds = readCDS(args.i + ".transdecoder.cds")
	writeBlocks(outfile, cds, index, args.p)

if __name__ == "__main__":
	main()
