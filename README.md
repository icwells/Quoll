# Quoll is a series of scripts for filtering and analyzing transcriptome alignments
Version 1.1 is meant to facilitate the analysis of transcriptomes assembled with 
Trinity and aligned with ProgressiveCactus (although Lastz/MultiZ alignments can 
also be used).

## Installation
Since most of the scripts are written in python3, simply download the repository:
git clone https://github.com/icwells/Quoll.git

The only script which needs to be compiled is mafTransToFasta.cpp. After downloading 
the repository, change into the new directory and enter the following:
g++ -o mafTransToFasta mafTransToFasta.cpp
(Depending on your system, you might use a C++ compiler other than g++).

## Getting Started
If you are starting with an alignment in maf format, first run mafTransToFasta which 
will extract DNA sequences from the maf file. Next, removeAncestor.py can be used to 
remove ancestral sequences from the alignment (if desired). Lastly, filterCDS.py can
be used to extract all transcripts with equivalent open reading frames from each 
aligned block. The ramaining scripts can then be called in any order for your analysis.

## Scripts
### mafTransToFasta
This script will pull DNA sequences from aligned blocks within a maf file and writes 
them in fasta format. Any sequences under 60 nucleotides in length will not be written 
since the shortest known protein is 20 amino acids long. The output file will be printed 
in the same directory as the input.

./mafTransToFasta {path to input maf}

### removeAncestor.py
This script will remove ancestral sequences which were inferred by alignment programs 
such as ProgressiveCactus.

python removeAncestor.py -i {path to input fasta alignment} -o {path to output file}

### filterCDS.py
This script identifies and writes common open reading frames within each aligned block. 
This is a fairly simple technique for extracting ORFs and will be replaced with a more 
sophisticated algorithm in the future.

python filterCDS.py -i {path to input fasta alignment} -o {path to output file} -l {minimum sequence length}

### alignedGenes.py
This script creates a csv table of equivalent gene IDs and any associated transcipt IDs 
for each alignment block. It also summarizes the total number of genes and transcripts 
found in each block.

python alignedGenes.py {path to input fasta alignment} 

### sampleMatrix.py
This script produces a csv matrix summarizing the number of genes aligned between each 
species/sample in the alignment. If "--combinations" is specified, it will return a matrix 
of the number of genes uniquely aligned between eavery possible combination of samples.

python sampleMatrix.py --combinations {path to input fasta alignment} 
 
### subsetAlignment.py
This script will subset aligned blocks if they have the specified number of samples aligned 
or if they contain genes from specified samples. The "--unique" flag indicates that blocks 
with only the given samples will be selected. Otherwise, blocks that contain at least those 
samples will be selected.

python subsetAlignment.py -n {# of samples to select} -i {path to input fasta alignment}
python subsetAlignment.py --unique -s {comma seperated list of samples to select} -i {path to input fasta alignment} 

### extractSequences.py
This script will extract one sequence from each aligned block. If no sample ID is given, 
the sequence will be chosen at random.

python extractSequences.py -t {sample to extract} -i {path to input fasta alignment} 

## Other Scripts
### transHeaders.py
This script contains functions for parsing gene and transcript IDs from transcriptomes 
and is called by several of the above scripts.

### cdsFinder.py
This script will identify all possible open reading frames for an alingment block. It 
is called by filterCDS.py but can also be called directly.

### blastToCSV.py
This script will combine blastn, blastx, and blastp results (in outfmt6 format) into a csv file. 
Input files must be titled blastn.outfmt6, blastx.outfmt6, and blastp.outfmt6.

python blastToCSV.py -i {path to directory containing blast results}  -o {path to output file}
