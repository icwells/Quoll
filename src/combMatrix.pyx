'''This script contains functions for calculating the number of unique
sample combinations in a trasncriptome alignment.'''

from collections import OrderedDict
from transHeaders import parseHeader, findSpecies

def extractIDs(ids):
	# Isolates sample ID from fasta header
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
	cdef str comb
	cdef str i
	# Assemble combination name
	samples.sort()
	comb = samples[0]
	for i in samples[1:]:
		comb += "-" + i
	if comb in counts.keys():
		# Increment if combination is present
		counts[comb] += 1
	else:
		# Add new entry
		counts[comb] = 1
	return counts		

def alignedSamples(infile):
	# Extracts sample IDs from each aligned block
	print("\tCalculating all combinations of aligned genes...\n")
	cdef counts = {}
	cdef list ids = []
	cdef str line
	with open(infile, "r") as alignment:
		for line in alignment:
			if line[0] == ">":
				ids.append(line[1:].strip())
			elif not line.strip():
				# Isolate sample IDS and count combinations
				samples = extractIDs(ids)
				if len(samples) >= 2:
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
	cdef str entry
	cdef list entries
	mat = OrderedDict()
	matrix = OrderedDict()
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
	# Sort data
	entries = list(mat.keys())
	entries.sort(key = len)
	for entry in entries:
		matrix[entry] = mat[entry]
	return matrix
