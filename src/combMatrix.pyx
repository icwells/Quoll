'''This script contains functions for calculating the number of unique
sample combinations in a trasncriptome alignment.'''

import itertools 
from collections import OrderedDict
from transHeaders import parseHeader, findSpecies

def makeComb(samples):
	# Calls itertools and converts result to list
	cdef list comb
	cdef list combs = [] 
	cdef str c
	cdef str j
	cdef int i
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
	cdef str i
	counts = OrderedDict()
	combs = makeComb(species)
	for i in combs:
		# Add entry
		counts[i] = 0
	return counts

def extractIDs(ids):
	# Isolates gsample ID from fasta header
	cdef list samples = []
	cdef str i
	for i in ids:
		name, delim = parseHeader(i, True)
		if name not in samples:
			# Add unique sample IDs to list
			samples.append(name)
	return samples

def countSamples(counts, samples):
	# Counts the number of matches for each category
	cdef str i
	combs = makeComb(samples)
	for i in combs:
		if i in counts.keys():
			# Increment if combination is present
			counts[i] += 1
	return counts		

def alignedSamples(infile, counts):
	# Extracts sample IDs from each aligned block
	print("\tCalculating all combinations of aligned genes...\n")
	cdef list ids
	cdef str line
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
	cdef int l
	cdef list row = []
	cdef int i
	cdef str count
	cdef str head
	cdef str tail
	mat = OrderedDict()
	l = len(species)
	# Create template list
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
