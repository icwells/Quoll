'''This script will subset a fasta alignment to include only genes with a 
given number of samples or with specific samples  present.'''

import argparse
import os
from transHeaders import parseHeader, findSpecies

def totalSpecies(species, ids):
	# Analyze genes
	names = []
	for i in ids:
		# Split gene ID from isoform number
		trans, delim = parseHeader(i, False)
		sp = trans[:trans.find(delim)]
		if sp not in names:
			names.append(sp)
	# Return number of species
	return len(names)			

def subByNumber(species, infile, outfile, n):
	# Returns aligned genes with one target species
	print("\tSubsetting aligned genes by number...")
	l = len(species)
	if l < n:
		# Quit if too many genes are specified
		print("\tPlease specify less than or equal to the total number of \
species/samples in the alignment.")
		quit()
	count = 0
	ids = []
	seqs = []
	with open(infile, "r") as alignment:
		with open(outfile, "w") as output:
			for line in alignment:
				if not line.strip():
					# Write output if the corrent number of seqs is present
					s = totalSpecies(species, ids)
					if s == n:
						for i in seqs:
							output.write(i)
						output.write("\n")
						count += 1
					# Clear lists
					ids = []
					seqs = []
				else:
					# Save formatted lines and gene IDs
					seqs.append(line)
					if line[0] == ">":
						ids.append(line[1:].strip())
	print(("\t{} genes with {} species aligned.\n").format(count, n))

#-----------------------------------------------------------------------------

def samplesPresent(species, ids, keep, unique):
	# Analyze genes
	names = []
	save = True
	for i in ids:
		# Split gene ID from isoform number
		trans, delim = parseHeader(i, False)
		sp = trans[:trans.find(delim)]
		if sp not in names:
			names.append(sp)
	# Return false if mismatches are found
	if keep == True:
		if unique == True:
			if len(names) != len(species):
				# Skip genes with additional samples
				return False
		else:
			if len(names) < len(species):
				# Skip genes if they are too short
				return False
		for i in species:
			if i not in names:
				# Don't print gene if any names are not present
				return False
	elif keep == False:
		for i in names:
			if i not in species:
				# Don't print gene if any names are present
				return False		
	# Return true if mismatches are not found
	return True

def subByName(infile, outfile, species, keep, unique, delim):
	# Returns aligned genes with one target species
	print("\tSubsetting aligned genes by sample ID...")
	count = 0
	ids = []
	seqs = []
	# Split species and sort
	species = species.split(",")
	with open(infile, "r") as alignment:
		with open(outfile, "w") as output:
			for line in alignment:
				if not line.strip():
					# Write output if the corrent number of seqs is present
					save = samplesPresent(species, ids, keep, unique)
					if save == True:
						for i in seqs:
							output.write(i)
						output.write("\n")
						count += 1
					# Clear lists
					ids = []
					seqs = []
				else:
					# Save formatted lines and gene IDs
					seqs.append(line)
					if line[0] == ">":
						ids.append(line[1:].strip())
	if keep == True:
		print(("\t{} genes with given samples aligned.\n").format(count))
	elif keep == False:
		print(("\t{} genes without given samples aligned.\n").format(count))

#-----------------------------------------------------------------------------

def main():
	parser = argparse.ArgumentParser(description = 
"This script will subset a fasta alignment to include only genes with a \
given number of samples or with specific samples  present.")
	parser.add_argument("-i", help = "Path to alignment file.")
	parser.add_argument("-n", type = int,
help = "Number of species required in alignment.")
	parser.add_argument("-s", help = "Comma seperated list of sample IDs to \
subset alignment by (exclusive of -n).")
	parser.add_argument("--unique", action = "store_true", 
help = "Output genes with only the given samples present (default is to \
output all genes with all samples present regardless of additional samples).")
	parser.add_argument("--exclude", action = "store_false",
help = "Indicates that the program should keep genes that do not include \
specified samples (default is to keep genes with specified samples).")
	args = parser.parse_args()
	species, delim = findSpecies(args.i)
	if args.n:
		# Make output file and subset
		outfile = args.i[:args.i.find(".")] + "-" + str(args.n) + "samples.fa"
		subByNumber(species, args.i, outfile, args.n)
	elif args.s:
		# Make output file and subset
		outfile = args.i[:args.i.find(".")] + "-sampleSubset.fa"
		subByName(args.i, outfile, args.s, args.exclude, args.unique, delim)

if __name__ == "__main__":
	main()
