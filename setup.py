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
		define_macros=[("USE_COEFFICIENTS", 1), ("PYTHON_EXTENSION", 1), ("NDEBUG", 1), ("ASSEMBLE_REDUCTION_MATRIX", 1)],
		extra_compile_args = options,
		language="c++",
        include_dirs=[numpy.get_include()]
	)),
	install_requires=[
        'matplotlib'
      ],
	name="ripser",
	description="Python wrapper around Uli Bauer's ripser code",
	author="Chris Tralie",
	author_email="chris.tralie@gmail.com",
	version='0.1',
	py_modules=['ripser'],
	# include_dirs=numpy.get_include()
)
