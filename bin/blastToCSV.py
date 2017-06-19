'''This script will concatenate blastn, blastx, and blastp results into 
a single csv summary annotation.'''

import argparse

def blastDict(indir):
	# Builds dictionary from blast files
	annotation = {}
	blastn = indir + "blastn.outfmt6"
	blastx = indir + "blastx.outfmt6"
	blastp = indir + "blastp.outfmt6"
	annotation = readBlast(annotation, blastp, "p")
	annotation = readBlast(annotation, blastx, "p")
	annotation = readBlast(annotation, blastn, "n")
	return annotation

def readBlast(annotation, blast, typ):
	# Append blast resutls to dict
	anno = {}
	with open(blast, "r") as bl:
		for line in bl:
			splt = line.split("\t")
			trans = splt[0]
			if "::" in trans:
				# Extract transcript from transdecoder IDs
				trans = trans.split("::")[1]
			name = splt[1]
			e = splt[10]
			if trans in anno.keys():
				if name == anno[trans][0]:
					# Skip if gene/protein names are the same
					pass
				else:
					# Replace if e-value is lower
					if e < anno[trans][1]:
						anno[trans] = [name, e]
			else:
				anno[trans] = [name, e]
	if typ == "p":
		# Determine correct protein name
		annotation = checkProtein(annotation, anno)
	elif typ == "n":
		for trans in annotation:
			if trans in anno.keys():
				# Add data to dict
				annotation[trans].extend([anno[trans][0], anno[trans][1]])
			else:
				# Add "NA" for missing values
				annotation[trans].extend(["NA", "NA"])
	return annotation

def checkProtein(annotation, anno):
	# Resolves conflicts between blastp and blastx proteins
	for trans in anno:
		if trans in annotation.keys():
			# Identify mismathces
			if annotation[trans][0] != anno[trans][0]:
				# Select protein with lower e-value
				if annotation[trans][1] > anno[trans][1]:
					annotation[trans][0] = anno[trans][0]
					annotation[trans][1] = anno[trans][1]
		else:
			# Append new protein entry
			annotation[trans] = anno[trans]
	return annotation

#-----------------------------------------------------------------------------

def writeCSV(annotation, outfile):
	# Prints dictionary to csv
	with open(outfile, "w") as output:
		output.write("TranscriptID,Protein,e,geneName,e\n")
		for line in annotation:
			# Add dict elements to string
			string = line + ","
			for i in annotation[line]:
				string += i + ","
			# Replace last comma with newline
			string = string[:-1] + "\n"
			output.write(string)

def main():
	parser = argparse.ArgumentParser(description = "This script will \
concatenate blastn, blastx, and blastp results into a single csv summary.")
	parser.add_argument("i", help = "Path to directory containing blast \
results. Must contain three files titled: blastn.outfmt6, blastx.outfmt6, \
and blastp.outfmt6. Output will be written to this directory.")
	args = parser.parse_args()
	indir = ars.i
	if indir[-1] != "/":
		indir += "/"
	outfile = indir + "mergedBlastResults.csv" 
	annotation = blastDict(args.i) 
	writeCSV(annotation, args.o)

if __name__ == "__main__":
	main()
