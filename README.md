# Quoll is a series of scripts for filtering and analyzing transcriptome alignments
Version 0.2 is meant to facilitate the analysis of transcriptomes assembled with 
Trinity and aligned with ProgressiveCactus (although Lastz/MultiZ alignments can 
also be used). 

Copyright 2017 by Shawn Rupp

## Installation
Download the repository:

git clone https://github.com/icwells/Quoll.git

Change into the Quoll directory:

cd Quoll/

Most of the scripts are written in python3, but several contain Cython modules which
must be compiled. Additionally, Quoll requires TransDecoder to identify open reading 
frames within transcripts. Running the install script will compile the Cython scripts 
and install TransDecoder from GitHub:

./install.sh

## Getting Started
If you are starting with an alignment in maf format, first run mafTransToFasta which 
will extract DNA sequences from the maf file.

cd bin/

python mafTransToFasta.py {path to input file}

Next, removeAncestor.py can be used to remove ancestral sequences from the alignment (if desired).

python removeAncestor.py {path to input file}

Lastly, filterCDS.py can be used to extract ORFs from the aligned transcripts. First, 
TransDecoder must be called on the output of the previous script to identify open reading 
frames. Export the path to TransDecoder and change into the output directory (TransDecoder 
writes output to the working directory).

export PATH=$PATH:{path to TransDecoder directory}

cd {output directory}

TransDecoder.LongOrfs -t {input file}

TransDecoder.LongOrfs -t {input file}

filterCDS.py can now be called to convert the TransDecoder output back into a fasta 
alignment. The input file for this script is the same input used for TransDecoder 
(its output will be inferred from the input file name). This script only needs the 
transdecoder.cds file, so all other TransDecoder output can be deleted if desired.

cd ~/Quoll/bin/

python filterCDS.py {path to input file}

The remaining scripts can then be called in any order for your analysis.
