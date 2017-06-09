'''This function will convert a maf transcriptome alignment to fasta format.'''

import os

def convert(infile, outfile):
	# Reads input and writes converted output to file
	cdef list block = []
	cdef str i
	cdef str line
	cdef list seq
	cdef str header
	cdef str nucl
	cdef str delim
	cdef int n = 1
	cdef int b = 1
	cdef int first = 0
	with open(infile, "r") as maf:
		with open(outfile, "w") as fasta:
			for line in maf:
				if not line.strip():
					# Write alinged block at first empty line
					if first == 1:
						for i in block:
							fasta.write(i)
						fasta.write("\n")
						b += 1
						first = 0
					block = []
				elif line[0] == "s":
					# Identify delimiter and isolate header and sequence
					delim = "\t"
					if " " in line:
						if line.find("\t") > line.find(" "):
							delim = " "
					seq = line.strip().split(delim)
					header = seq[1].strip()
					# Append block and gene index number to header
					header += ("|{}::{}").format(b, n)
					nucl = seq[-1].strip()
					if len(nucl) >= 60:
						# Write header and sequence to fasta
						block.append((">{}\n{}\n").format(header, nucl))
						n += 1
						first = 1
