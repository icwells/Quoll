'''This script will provide a list of genes which have aligned between 
different species.'''

import argparse
import os
from transHeaders import parseHeader, findSpecies

def assemble(species, genes):
	# Converts dict to ordered string
	string = ""
	for i in species:
		try:
			# Add entries for aligned species
			for j in genes[i]:
				string += j + "; "
			string = string[:-2] + ","
		except KeyError:
			# Add dash for missing species
			string += "-,"
	return string

def findOrthologs(species, ids, idtype):
	# Analyze genes
	genes = {}
	t = len(ids)
	for i in ids:
		# Split gene ID from isoform number
		trans, delim = parseHeader(i, False)
		sp = trans[:trans.find(delim)]
		iso = ""
		if trans.count(delim) > 1:
			geneid = trans[:trans.rfind(delim)]
			iso = trans[trans.rfind(delim)+1:]
		else:
			# Store refseq ID with no isoforms
			geneid = trans[trans.rfind(delim)+1:]
		if sp in genes.keys():
			new = True
			for idx,j in enumerate(genes[sp]):
				if geneid in genes[sp][idx]:
					genes[sp][idx] = genes[sp][idx][:-1] + " " + iso + ")"
					new = False
					break
			if new == True:
				genes[sp].append(geneid + " (" + iso + ")")
		else:
			# Save new entry as list item
			if iso:
				genes[sp] = [geneid + " (" + iso + ")"]
			else:
				# Store refseq ID with no isoforms
				genes[sp] = [geneid]
	# Build output string
	s = len(genes)
	ortho = assemble(species, genes)
	ortho += ("{},{}\n").format(t, s)
	return ortho, s					

def comparison(species, infile, outfile, idtype):
	# Returns aligned genes with one target species
	print("\tAnalyzing aligned genes...")
	if idtype == "tr":
		delim = "-"
	l = len(species)
	# Make list for number of transcripts + number of aligned genes 
	stats = [0]
	for i in range(0, l+1):
		stats.append(0)	
	n = 0
	ids = []
	with open(infile, "r") as alignment:
		with open(outfile, "w") as output:
			# Write header
			head = ""
			for i in species:
				head += i + ","
			output.write(head + "totalTranscripts,alignedGenes\n")
			for line in alignment:
				if not line.strip():
					# Write output line and collect stats
					ortho, s = findOrthologs(species, ids, idtype)
					output.write(ortho)
					# Append 1 to # of genes and column of aligned species 
					stats[0] += 1
					stats[s] += 1
					# Reset gene values
					n = 0
					ids = []
				elif line[0] == ">":
					ids.append(line[1:].strip())
					n += 1
	return stats

def printStat(species, stats):
	# Prints numbers of aligned genes
	l = len(species)
	print("\n\tSummary of alignment data:")
	print(("\tNumber of samples identified: {}").format(l))
	# Calls comparison function and prints species stats
	print(("\tNumber of genes aligned: {}\n").format(stats[0]))
	for idx,n in enumerate(stats[2:-1]):
		# Append number of genes present in n species 
		string = ("{} samples: {}").format(idx+2, n)
		print("\tNumber of genes present in " + string)

def main():
	parser = argparse.ArgumentParser(description = "This script will provide \
a list of genes which have aligned between different species.")
	parser.add_argument("i", help = "Path to alignment file.")
	args = parser.parse_args()
	# Make output directory
	outfile = args.i[:args.i.find(".")] + "-alignedGenes.csv"
	species, idtype = findSpecies(args.i)
	stats = comparison(species, args.i, outfile, idtype)
	printStat(species, stats)

if __name__ == "__main__":
	main()
