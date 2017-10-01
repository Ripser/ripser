from distutils.core import setup, Extension
import numpy

c_ext = Extension("ripser", ["_ripser.cpp", "ripser.cpp"])

setup(
    ext_modules=[c_ext],
    include_dirs=numpy.get_include(),
    define_macros=[("USE_COEFFICIENTS", 1), ("ASSEMBLE_REDUCTION_MATRIX", 1), ("PYTHON_EXTENSION", 1)],
    extra_compile_args = ["-std=c++11"]
)
