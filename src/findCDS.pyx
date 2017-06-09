'''These functions will sort TransDecoder output to identify open reading frames 
in a nucleotide alignment and remove genes with insufficient coverage.'''

from collections import OrderedDict

def writeBlocks(outfile, cds, index, percent):
	# Writes reformated TransDecoder output in fasta alignment format
	print("\tWriting cds alignment...")
	cdef int b = 0
	cdef int n = 0
	cdef double content
	cdef str block
	cdef str i
	with open(outfile, "w") as output:
		for block in index:
			if block in cds.keys():
				for i in index[block]:
					if i in cds[block].keys():
						content = 1 - (cds[block][i][1].count("-")/
										len(cds[block][i][1].strip()))
						if content >= percent:
							# Only write if nucleotide content is above min
							n += 1
							output.write(cds[block][i][0] + cds[block][i][1] + "\n")
				# Add newline between blocks
				b += 1
				output.write("\n")
	print(("\tWrote {:,} blocks and {:,} transcripts.\n").format(b, n))

def trimHeader(line):
	# Trims excess info from TransDecoder headers
	cdef str header
	cdef list block
	line = line.split()[1]
	line = line[line.find("::")+2:line.rfind("::")]
	line = line.split("|")
	header = (">{}\n").format(line[0])
	block = line[1].split("::")
	# Return gene ID and gene #
	return header, block[0], block[1] 

def readCDS(infile):
	# Edits TransDecoder output
	print("\tReading TransDecoder output...")
	cds = {}
	cdef int save = 0
	cdef int n = 0
	cdef str line
	cdef str seq
	cdef str gene
	with open(infile, "r") as fasta:
		for line in fasta:
			if line[0] == ">":
				if save == True:
					if block not in cds.keys():
						# Add entry for new block
						cds[block] = {}
					if gene in cds[block].keys():
						#print("\tMultiples of " + name)
						if len(seq) > len(cds[block][gene][1]):
							cds[block][gene] = [header, seq]
					else:
						cds[block][gene] = [header, seq]
				header, block, gene = trimHeader(line)
				seq = ""
			else:
				seq += line.strip()
				save = True
	return cds

def indexAlignment(infile):
	# Record which genes are aligned in each block
	print("\n\tIndexing fasta alignment...")
	index = OrderedDict()
	block = OrderedDict()
	cdef int n = 0
	cdef int b = 0
	cdef str line
	cdef list header
	cdef list idx
	with open(infile, "r") as fasta:
		for line in fasta:
			if line[0] == ">":
				# Isolate index numbers and append to dict
				line = line.strip()
				header = line.strip().split("|")
				idx = header[1].split("::")				
				if idx[0] not in index.keys():
					# Add new block
					index[idx[0]] = {}
				block[idx[1]] = header				
			elif not line.strip():
				index[idx[0]] = block
				block = OrderedDict()
	return index
