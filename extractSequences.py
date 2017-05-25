'''This script will will extract sequences from an alignment.'''

import argparse
from random import choice
from transHeaders import parseHeader

def extract(seqs, target):
	# Pull target from aligned black
	seq = {}		
	if target:
		for header in seqs:
			sample = parseHeader(header, True)
			if sample ==  target:
				# Save any matching entry
				seq[header] = seqs[header]
	else:
		# Extract random entry
		keys = list(seqs.keys())
		sample = choice(keys)
		seq[sample] = seqs[sample]
	return seq		

def extractFromFasta(infile, outfile, target):
	# Extracts given/random samples from each alignemnt block
	seqs = {}
	with open(infile, "r") as alignment:
		with open(outfile, "w") as output:
			for line in alignment:
				if not line.strip():
					# Write output if the corrent number of seqs is present
					seq = extract(seqs, target)
					if seq:
						for sample in seq:
							output.write((">{}\n{}\n").format(sample, 
														seq[sample]))
					# Clear dict
					seqs = {}
				elif line[0] == ">":
					header = line[1:].strip()
				else:
					# Save formatted lines and gene IDs
					seqs[header] = line

def main():
	parser = argparse.ArgumentParser(description = 
"This script will will extract sequences from an alignment. Outputs \
sequences in standard multi-fasta format (i.e. not an alignment)")
	parser.add_argument("-i", help = "Path to input fasta alignment.")
	parser.add_argument("-t", 
help = "Sample ID of target sequence to extract. If none is given, one \
sample will be randomly selected for each aligned block.")
	args = parser.parse_args()
	outfile = args.i[:args.i.find(".")] + "-ExtractedSequences.fa"
	extractFromFasta(args.i, outfile, args.t)

if __name__ == "__main__":
	main()
