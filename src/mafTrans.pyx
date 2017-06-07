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
	with open(infile, "r") as maf:
		with open(outfile + ".temp", "w") as fasta:
			for line in maf:
				if not line.strip():
					# Write alinged block
					for i in block:
						fasta.write(i)
					fasta.write("\n")
					block = []
				elif line[0] == "s":
					# Identify delimiter and isolate header and sequence
					delim = "\t"
					if " " in line:
						if line.find("\t") > line.find(" "):
							delim = " "
					seq = line.strip().split(delim)
					header = seq[1].strip()
					nucl = seq[-1].strip()
					if len(nucl) >= 60:
						# Write header and sequence to fasta
						block.append((">{}\n{}\n").format(header, nucl))
	rmEmpty(outfile)

def rmEmpty(outfile):
	# Removes consecutive empty lines in fasta alignment
	cdef int first = 0
	cdef str line
	with open(outfile + ".temp", "r") as temp:
		with open(outfile, "w") as output:
			for line in temp:
				if not line.strip():
					if first == 1:
						output.write("\n")
						first = 0
				else:
					output.write(line)
					first = 1
	os.remove(outfile + ".temp")

