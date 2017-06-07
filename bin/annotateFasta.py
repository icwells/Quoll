'''This script will add gene IDs and protein and gene names to fasta 
headers using blastToCSV output.'''

import argparse

def readDict(anno):
	# Converts blastToCSV output to dictionary of headers
	headers = {}
	with open(anno, "r") as annotation:
		for line in annotation:
			splt = line.split(",")
			# Gene ID, trnascript ID, protein name, and gene name
			headers[splt[0]] = (">" + splt[1] + "-" + splt[0] + "-" + 
								splt[2] + "-" + splt[4] + "\n")
	return headers

def annotate(infile, headers):
	# Replaces old header with annotated headers
	# Create output file
	if ".fasta" in infile:
		outfile = infile.replace(".fasta", ".annotated.fasta")
	else:
		outfile = infile.replace(".fa", ".annotated.fa")
	with open(infile, "r") as fasta:
		with open(outfile, "w") as output:
			for line in fasta:
				if line[0] == ">":
					trans = line.replace(">", "").strip()
					if trans in headers.keys():
					# Replace header if there is a match
						output.write(headers[trans])
					else:
						# Write original if there is no match
						output.write(line)
				else:
					# Write sequence
					output.write(line)					

def main():
	parser = argparse.ArgumentParser(description = "This script will add gene \
IDs and protein and gene names to fasta headers using blastToCSV output..")
	parser.add_argument("-i", help = "Path to input fasta")
	parser.add_argument("-a", help = "Path to blastToCSV annotation.")
	args = parser.parse_args()
	headers = readDict(args.a) 
	annotate(args.i, headers)

if __name__ == "__main__":
	main()
