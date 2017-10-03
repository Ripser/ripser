from distutils.core import setup, Extension
import numpy
import sys

#c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D PRINT_PERSISTENCE_PAIRS
options = ["-std=c++11", '-Ofast']

if sys.version_info[0] == 2:
	options.append("-fpermissive")

c_ext = Extension("_ripser",
sources = ["ripser.cpp", "_ripser.cpp"],
define_macros=[("USE_COEFFICIENTS", 1), ("ASSEMBLE_REDUCTION_MATRIX", 1), ("PYTHON_EXTENSION", 1), ("NDEBUG", 1)],
extra_compile_args = options,
)

setup(name="ripser",
    ext_modules=[c_ext],
    include_dirs=numpy.get_include()
)
