'''This script will remove ancestor sequences from a multiple species 
alignment.'''

import argparse

def rmAncestor(infile, outfile):
	rm = False
	newline = True
	with open(infile, "r") as fasta:
		with open(outfile, "w") as output:
			for line in fasta:
				if line == "\n":
					if newline == True:
						# Only write one blank line between genes
						output.write("\n")
						newline = False
				elif line[0] == ">":
					if line[:4] == ">Anc":
						# Skip ancestor headers
						rm = True
					else:
						# Write out species headers
						line = line.replace("_","-")
						output.write(line)
						newline = True
				else:
					if rm == False:
						# Write out species sequences
						output.write(line)
					elif rm == True:
						# Skip ancestor sequences
						rm = False		

def main():
	parser = argparse.ArgumentParser(description="This script will remove \
ancestor sequences from a multiple species alignment.")
	parser.add_argument("-i", help="Path to input file.")
	parser.add_argument("-o", help="Path to output file.")
	args = parser.parse_args()
	rmAncestor(args.i, args.o)

if __name__ == "__main__":
	main()
