'''This script will will extract aligned blocks from an alignment.'''

import argparse
from collections import OrderedDict
from random import choice
from transHeaders import parseHeader

def extract(seqs, targets):
	# Pull target from aligned black
	seq = OrderedDict()
	found = []
	for header in seqs:
		sample, delim = parseHeader(header, True)
		if sample in targets and sample not in found:
			# Save one matching entry for each sample
			seq[header] = seqs[header]
			found.append(sample)
	for i in targets:
		if i not in found:
			# Clear dict if sample is missing
			seq = {}
			break
	return seq		

def extractFromFasta(infile, outfile, targets):
	# Extracts given/random samples from each alignemnt block
	seqs = OrderedDict()
	with open(infile, "r") as alignment:
		with open(outfile, "w") as output:
			for line in alignment:
				if not line.strip():
					# Write output if the corrent number of seqs is present
					seq = extract(seqs, targets)
					if seq:
						for sample in seq:
							output.write((">{}\n{}").format(sample, 
														seq[sample]))
						output.write("\n")
					# Clear dict
					seqs = OrderedDict()
				elif line[0] == ">":
					header = line[1:].strip()
				else:
					# Save formatted lines and gene IDs
					seqs[header] = line

def main():
	parser = argparse.ArgumentParser(description = 
"This script will will extract  aligned blocks from an alignment. Only one \
sequence per species will be saved.")
	parser.add_argument("-i", help = "Path to input fasta alignment.")
	parser.add_argument("-t", 
help = "Comma-seperated list of sample IDs of target sequences to extract.")
	args = parser.parse_args()
	outfile = args.i[:args.i.find(".")] + "-sampleAlignment.fa"
	targets = args.t.split(",")
	extractFromFasta(args.i, outfile, targets)

if __name__ == "__main__":
	main()
