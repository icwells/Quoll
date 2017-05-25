'''This script will create a matrix of the number of genes aligned between 
every sample in an alignment.'''

import argparse
import os
import itertools
from collections import OrderedDict
from transHeaders import parseHeader, findSpecies

#-----------------------------------------------------------------------------

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

def makeComb(samples):
	# Calls itertools and converts result to list
	comb = []
	combs = []
	for i in range(len(samples)):
		if i > 1:
			comb = itertools.combinations(samples, i)
			for item in comb:
				c = ""
				for j in item:
					# Join sample names for combination name
					c += str(j) + "-"
				combs.append(c[:-1])
	return combs	

def combDict(species):
	# Creates entries for each combination
	counts = OrderedDict()
	combs = makeComb(species)
	for i in combs:
		# Add entry
		counts[i] = 0
	return counts

def extractIDs(ids):
	# Isolates gsample ID from fasta header
	samples = []
	for i in ids:
		name, delim = parseHeader(i, True)
		if name not in samples:
			# Add unique sample IDs to list
			samples.append(name)
	return samples

def countSamples(counts, samples):
	# Counts the number of matches for each category
	combs = makeComb(samples)
	for i in combs:
		if i in counts.keys():
			# Increment if combination is present
			counts[i] += 1
	return counts		

def alignedSamples(infile, counts):
	# Extracts sample IDs from each aligned block
	print("\tCalculating all combinations of aligned genes...\n")
	ids = []
	with open(infile, "r") as alignment:
		for line in alignment:
			if line[0] == ">":
				ids.append(line[1:].strip())
			elif not line.strip():
				# Isolate sample IDS and count combinations
				samples = extractIDs(ids)
				counts = countSamples(counts, samples)
				ids = []
	return counts

def convertToMatrix(counts, species):
	# Convert dict to printable matrix
	mat = OrderedDict()
	l = len(species)
	# Create template list
	row = []
	for i in range(l):
		row.append(0)
	for count in counts:
		# Split first sample from list and index species for column
		head = count[:count.find("-")]
		tail = count[count.find("-")+1:]
		idx = species.index(head)
		if tail in mat.keys():
			mat[tail][idx] = counts[count]
		else:
			# Save data in list entry
			mat[tail] = []
			mat[tail].extend(row)
			mat[tail][idx] = counts[count]
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
