from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy 

sourcefiles  = ['utils.pyx','sequenceAlign.cpp']
compile_opts = ['-std=c++14',"-O3"]
libs = ["gsl", "gslcblas"]
includes = [".", "/scicore/home/neher/GROUP/bin/anaconda2/envs/cpp/include",numpy.get_include()]

ext=[Extension('*',
            sources=sourcefiles,
            libraries=libs,
            extra_compile_args=compile_opts,
            include_dirs=includes,
            language='c++')]

setup(
  ext_modules=cythonize(ext)
  )