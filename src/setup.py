'''This script will cythonize the Quoll package.'''

from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("combMatrix.pyx"))
setup(ext_modules=cythonize("findCDS.pyx"))
setup(ext_modules=cythonize("mafTrans.pyx"))
setup(ext_modules=cythonize("transHeaders.pyx"))
