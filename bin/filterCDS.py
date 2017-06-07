'''This program will sort TransDecoder output to identify open reading frames 
in a nucleotide alignment and remove genes with insufficient coverage.'''

import shutil
import os
import argparse

def writeBlocks(outfile, cds, index, percent):
	# Writes reformated TransDecoder output in fasta alignment format
	print("\tWriting cds alignment...")
	b = 0
	n = 0
	with open(outfile, "w") as output:
		for block in index:
			for i in block:
				if i in cds.keys():
					content = 1 - (cds[i][1].count("-")/len(cds[i][1].strip()))
					if content >= percent:
						# Only write if nucleotide content is above min
						n += 1
						output.write(cds[i][0] + cds[i][1] + "\n")
			# Add newline between blocks
			b += 1
			output.write("\n")
	print(("\tWrote {:,} blocks and {:,} transcripts.").format(b, n))

def trimHeader(line):
	# Trims excess info from TransDecoder headers
	line = line.split()[1]
	return line.split("::")[1]

def readCDS(infile):
	# Edits TransDecoder output
	print("\tReading TransDecoder output...")
	cds = {}
	save = False
	n = 0
	with open(infile, "r") as fasta:
		for line in fasta:
			if line[0] == ">":
				if save == True:
					n += 1
					if name in cds.keys():
						#print("\tMultiples of " + name)
						if len(seq) > len(cds[name][1]):
							cds[name] = [header, seq]
					else:
						cds[name] = [header, seq]
				name = trimHeader(line)
				header = (">{}\n").format(name)
				seq = ""
			else:
				seq += line.strip()
				save = True
	print(("\tIdentified {:,} open reading frames.").format(n))
	return cds

def indexAlignment(infile):
	# Record which genes are aligned in each block
	print("\n\tIndexing fasta alignment...")
	index = []
	block = []
	n = 0
	b = 0
	with open(infile, "r") as fasta:
		for line in fasta:
			if not line.strip():
				# Save IDs from each block
				index.append(block)
				block = []
				b += 1
			elif line[0] == ">":
				block.append(line[1:].strip())
				n += 1
	print(("\tIdentified {:,} blocks and {:,} transcripts.").format(b, n))
	return index	
			 
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
