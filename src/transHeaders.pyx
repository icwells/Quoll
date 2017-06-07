'''This script contains redundant functions for transcriptome headers.'''

import re

def trinityIDs(header, idonly):
	# Returns species ID or shortened geneid from Trinity header line
	cdef str name
	cdef str delim
	header = header.replace(">", "").strip()
	name = header.split("TRINITY")[1][1:]
	delim = "_"
	if "-" in name:
		delim = "-"
	if idonly == True:
		name = name.split(delim)[0]
	return name, delim

def refseqIDs(header, idonly):
	# Returns species ID or geneid from NCBI refSeq header line
	cdef str name
	name = header.replace(">", "")
	if idonly == True:
		name = name.split(".")[0]
	return name, "."

def parseHeader(header, idonly):
	# Sends header to appropriate parsing function (only returns name)
	cdef str name
	cdef str delim
	if "TRINITY" in header:
		name, delim = trinityIDs(header, idonly)
	elif re.search(r".rna(\d+)", header):
		name, delim = refseqIDs(header, idonly)
	return name, delim

def findSpecies(infile):
	# Gathers species names from input
	cdef list species = []
	cdef str line
	cdef str name
	cdef str delim
	print("\n\tIdentifying species in alignment...")
	with open(infile, "r") as alignment:
		for line in alignment:
			if line[0] == ">":
				if "TRINITY" in line:
					name, delim = trinityIDs(line, True)
				elif re.search(r".rna(\d+)", line):
					name, delim = refseqIDs(line, True)
				if name not in species:
					species.append(name)
	species.sort()
	return species, delim
