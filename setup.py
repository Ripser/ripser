from distutils.core import setup, Extension
import numpy
import sys

options = ["-std=c++11"]

if sys.version_info[0] == 2:
	options.append("-fpermissive")

c_ext = Extension("ripser", 
sources = ["ripser.cpp", "_ripser.cpp"],
define_macros=[("USE_COEFFICIENTS", 1), ("ASSEMBLE_REDUCTION_MATRIX", 1), ("PYTHON_EXTENSION", 1)],
extra_compile_args = options,
)

setup(name="ripser",
    ext_modules=[c_ext],
    include_dirs=numpy.get_include()
)
