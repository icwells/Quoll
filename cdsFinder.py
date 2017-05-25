'''This script will identify putative open reading frames in a single gene 
alignment and write the statistics to an output file.'''

from collections import OrderedDict
import argparse

def seqDict(infile):
	# Removes seqs with low content, identifies starts and stops
	seqs = OrderedDict()
	with open(infile, "r") as fasta:
		for line in fasta:
			if line[0] == ">":
				species = line[1:].strip()
			else:
				seqs[species] = line.upper().strip()
	return seqs

def cdsStats(seqs):
	# Only proceed if nucleotides compose greater than the percent
	stats = OrderedDict()
	for species in seqs:
		# Get count of all nucleotides
		seq = seqs[species]
		l = len(seq)
		aligned = (l - seq.count("-"))/l*100
		stats[species] = [l, aligned]
		idx = 0
		while idx < 3:
			# Examine each ORF
			orf = startStop(seq, idx)
			stats[species].append(orf)
			idx += 1
	return stats

def startStop(seq, idx):
	# Identify all possible start and stop codons
	codons = []
	starts = []
	stops = []
	seq = seq[idx:]
	for i in range(0, len(seq), 3):
		try:
			codons.append(seq[i:i +3])
			i += 3
		except IndexError:
			# Skip trailing nucleoltides
			pass
	# Identify in frame starts
	for ind,i in enumerate(codons):
		if i == "ATG":
			starts.append(ind*3)
		elif i == "TAA" or i == "TAG" or i == "TGA":
			stops.append(ind*3)
	# Append empty value if starts/stops not found
	if not starts:
		starts = ["NA"]
	if not stops:
		stops = ["NA"]
	return [starts, stops]

def writeDict(stats, outfile):
	# Write stats sequences to file
	with open(outfile, "w") as output:
		output.write("Possible start and stop codons for each sample\n\n")
		for species in stats:
			idx = 1
			output.write(species + "\n")
			output.write(("\tTranscript length: {}\n").format(stats[species][0]))
			output.write(("\tPercent of aligned nucleotides: {}\n").format(
							stats[species][1]))
			output.write("\tPossible start and stop codon indecies:\n")
			while idx < 4:
				output.write(("\t\tORF {}:\n\t\t\tStart Codons:").format(idx))
				for i in stats[species][idx+1][0]:
					output.write(" " + str(i))
				output.write("\n\t\t\tStop Codons:")
				for i in stats[species][idx+1][1]:
					output.write(" " + str(i))
				output.write("\n")
				idx += 1
			output.write("\n")
			 
def main():
	parser = argparse.ArgumentParser(description = 
"This script will identify putative open reading frames in a single gene \
alignment and write the statistics to an output file")
	parser.add_argument("i", help = "Path to input file.")
	args = parser.parse_args()
	# Make output file
	outfile = args.i[:args.i.find(".")] + "-possibleORFs.txt"
	# Filter for content first to remove empty lines
	seqs = seqDict(args.i)
	stats = cdsStats(seqs)
	writeDict(stats, outfile)

if __name__ == "__main__":
	main()
