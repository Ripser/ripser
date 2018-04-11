import sys
from distutils.core import setup, Extension

from Cython.Build import cythonize
from Cython.Distutils import build_ext


#c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D PRINT_PERSISTENCE_PAIRS
#options = ["-std=c++11", "-Ofast"]
options = ["-std=c++11", "-Ofast"]

if sys.version_info[0] == 2:
	options.append("-fpermissive")

class CustomBuildExtCommand(build_ext):
    """build_ext command for use when numpy headers are needed."""
    def run(self):

        # Import numpy here, only when headers are needed
        import numpy

        # Add numpy headers to include_dirs
        self.include_dirs.append(numpy.get_include())

        # Call original build_ext command
        build_ext.run(self)

setup(ext_modules = cythonize(
	Extension("pyRipser",
		sources = ["src/pyRipser.pyx"],
		define_macros=[("USE_COEFFICIENTS", 1), ("PYTHON_EXTENSION", 1), ("NDEBUG", 1), ("ASSEMBLE_REDUCTION_MATRIX", 1)],
		extra_compile_args = options,
		language="c++"
	)),
	install_requires=[
		'Cython',
		'numpy',
		'scipy',
        'matplotlib',
		'scikit-learn'
    ],
	cmdclass = {'build_ext': CustomBuildExtCommand},
	name="ripser",
	description="Persistent homology for people",
	long_description="A Python wrapper of Ripser meant to integrate into your workflow.",
	author="Chris Tralie, Nathaniel Saul",
	author_email="chris.tralie@gmail.com, nathaniel.saul@wsu.edu",
	version='0.1.1',
	package_dir = {'': 'src'},
	py_modules=['ripser']
)
