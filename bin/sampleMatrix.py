'''This script will create a matrix of the number of genes aligned between 
every sample in an alignment.'''

import argparse
from collections import OrderedDict
from transHeaders import parseHeader, findSpecies
from combMatrix import combDict, alignedSamples, convertToMatrix

def matrix(species):
	# Create matrix 
	mat = OrderedDict()
	row = []
	for x in range(len(species)):
		row.append(0)
	for i in species:
		# Extend list to avoid setting entries equal to same variable
		mat[i] = []
		mat[i].extend(row)
	return mat		

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
	return names

def matrixAppend(mat, species, names):
	# Increments matrix entries	
	for name in names:
		for i in names:
			mat[name][species.index(i)] += 1
	return mat

def countAlignment(species, infile, delim):
	# Returns aligned genes with one target species
	print("\tCalculating number aligned genes...\n")
	mat = matrix(species)
	ids = []
	with open(infile, "r") as alignment:
		for line in alignment:
			if not line.strip():
				# Write output if the corrent number of seqs is present
				names = totalSpecies(species, ids)
				mat = matrixAppend(mat, species, names)
				# Clear lists
				ids = []
			elif line[0] == ">":
				ids.append(line[1:].strip())
	return mat

#-----------------------------------------------------------------------------

def printMatrix(mat, species, outfile):
	# Saves matrix as csv
	with open(outfile, "w") as output:
		row = ""
		for i in species:
			row += "," + i
		output.write(row + "\n")
		for name in mat:
			row = name
			for i in mat[name]:
				if i == 0:
					row += ",-"
				else:
					row += "," + str(i)
			output.write(row + "\n")

def main():
	parser = argparse.ArgumentParser(description = 
"This script will subset a fasta alignment to include only genes with a given \
number of species/samples present.")
	parser.add_argument("--combinations", action = "store_true",
help = "Returns matrix of genes aligning between all possible combinations of \
samples (By default it will return values for pairwise combinations)." )
	parser.add_argument("i", help = "Path to alignment file.")
	args = parser.parse_args()
	species, delim = findSpecies(args.i)
	if args.combinations == True:
		# Make output file
		outfile = args.i[:args.i.find(".")] + "-combinationMatrix.csv"
		counts = combDict(species)
		counts = alignedSamples(args.i, counts)
		mat = convertToMatrix(counts, species)
		printMatrix(mat, species, outfile)
	else:
		# Make output file
		outfile = args.i[:args.i.find(".")] + "-sampleMatrix.csv"
		mat = countAlignment(species, args.i, delim)
		printMatrix(mat, species, outfile)

if __name__ == "__main__":
	main()
