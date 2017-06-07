#!/bin/bash

##############################################################################
# This bash script will cythinize the Quoll package and install TransDecoder
# 	Note: By running this file, you are agreeing to the licensing terms of both
#		the Quoll package and TransDecoder. 
##############################################################################

cd src/
python setup.py build_ext --inplace
cd ../

mv src/combMatrix.*.so bin/combMatrix.so
mv src/mafTrans.*.so bin/mafTrans.so
mv src/transHeaders.*.so bin/transHeaders.so

git clone https://github.com/TransDecoder/TransDecoder.git
cd TransDecoder
make
