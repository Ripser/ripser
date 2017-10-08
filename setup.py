from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
import sys

#c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D PRINT_PERSISTENCE_PAIRS
#options = ["-std=c++11", "-Ofast"]
options = ["-std=c++11", "-Ofast"]

if sys.version_info[0] == 2:
	options.append("-fpermissive")

setup(ext_modules = cythonize(
	Extension("pyRipser",
		sources = ["pyRipser.pyx"],
		define_macros=[("USE_COEFFICIENTS", 1), ("PYTHON_EXTENSION", 1), ("NDEBUG", 1)],
		extra_compile_args = options,
		language="c++",
	)),
	include_dirs=numpy.get_include()
)
