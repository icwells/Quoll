'''This program will identify open reading frames in a nucleotide alignment 
and remove genes with insufficient coverage.'''

from cdsFinder import cdsStats, startStop
from transHeaders import parseHeader
from collections import OrderedDict
import argparse

def extractCDS(infile, outfile, percent, minlen, retain):
	# Removes seqs with low content, identifies starts and stops
	seqs = OrderedDict()
	with open(outfile, "w") as output:
		with open(infile, "r") as fasta:
			for line in fasta:
				if not line.strip():
					# Filter for content first to remove empty lines
					cont, seqs, orfs = findORFs(seqs, minlen, percent)
					if cont == True:
						wr, seqs = findConsensus(seqs, orfs, retain)
						if wr == True:
							writeDict(seqs, output)
					seqs = OrderedDict()
				elif line[0] == ">":
					species = line[1:].strip()							
				else:
					seqs[species] = line.upper().strip()
		# Write last alignment
		cont, seqs, orfs = findORFs(seqs, minlen, percent)
		if cont == True:
			seqs = findConsensus(seqs, orfs)
			writeDict(seqs, output)

#-----------------------------------------------------------------------------

def examineORFs(orfs, minlen):
	# Determine most likely ORF
	indecies = []
	for frame in orfs:	
		starts = frame[0]
		stops = frame[1]
		for i in starts:
			if i != "NA":
				for j in stops:
					if j != "NA":
						# Find length (add 3: stop is index of 1st nucleotide)
						l = (j - i) + 3
						if l % 3 == 0 and l >= minlen:
							# Append all possible ORFs
							indecies.append([l, i, j])
	return indecies

def findORFs(seqs, minlen, percent):
	# Calls cdsFinder and determine common ORF
	orfs = {}
	stats = cdsStats(seqs)
	for sample in stats:
		# Collect stats for each sample
		if stats[sample][1] < percent:
			# Skip transcripts with low content
			del seqs[sample]
		else:
			species, delim = parseHeader(sample, True)
			orf = examineORFs(stats[sample][2:], minlen)
			if orf:
				try:
					orfs[species][sample] = orf
				except KeyError:
					orfs[species] = {}
					orfs[species][sample] = orf
			else:
				del seqs[sample]
	if len(seqs.keys()) > 1:
		return True, seqs, orfs
	else:
		return False, seqs, orfs

def mode(lengths):
	# Determines mode(s) of list. Returns list
	counts = {}
	matches = []
	# Count total occurances
	for i in lengths:
		try:
			counts[i] += 1
		except KeyError:
			counts[i] = 0
			counts[i] += 1
	# Determine most common
	maximum = 0
	for i in counts:
		if counts[i] > maximum:
			maximum = i
			matches.append(i)
		elif counts[i] == maximum:
			matches.append(i)
	return matches

def findConsensus(seqs, orfs, retain):
	# Finds most common lengths
	lengths = []
	for species in orfs:
		lens = []
		for sample in orfs[species]:
			for i in orfs[species][sample]:
				lens.append(i[0])
		# Append species' most common length to avoid biasing overall estimate
		matches = mode(lens)
		lengths.extend(matches)
	# Find overall most common length
	l = mode(lengths)
	l = l[0]
	for species in orfs:
		for sample in orfs[species]:
			match = False
			for i in orfs[species][sample]:
				if i[0] == l:
					# Trim transcript (add 2 to stop for whole codon)
					seqs[sample] = seqs[sample][i[1]:i[2]+3]
					match = True
					break
			if match == False:
				if retain == True:
					# Keep l bases after first start
					idx = seqs[sample].find("ATG")
					seqs[sample] = seqs[sample][idx:idx+l]
				elif retain == False:
					del seqs[sample]
	if len(seqs.keys()) > 1:
		return True, seqs
	else:
		return False, seqs

#-----------------------------------------------------------------------------

def writeDict(seqs, output):
	# Write passing sequences to file
	for species in seqs:
		output.write(">" + species + "\n" + seqs[species] + "\n")
	output.write("\n")
			 
def main():
	parser = argparse.ArgumentParser(description = "This program will \
identify open reading frames in a nucleotide alignment and remove genes with \
insufficient coverage.")
	parser.add_argument("-i", help = "Path to input file.")
	parser.add_argument("-o", help = "Path to output file.")
	parser.add_argument("-p", type=float, default=0.9, help="Minimum required \
percentage of nucleotides remaining after filtering (as a decimal).")
	parser.add_argument("-l", type = int, default = 60,
 help = "Minimum transcript length (default = 60bp - smallest known protein).")
	parser.add_argument("--retain", action = "store_true", 
help = "Retains aligned seqeunces which do not have an ORF of matching length.")
	args = parser.parse_args()
	extractCDS(args.i, args.o, args.p, args.l, args.retain)

if __name__ == "__main__":
	main()
