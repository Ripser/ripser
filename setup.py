import sys

from setuptools import setup
from setuptools.extension import Extension

# we'd better have Cython installed, or it's a no-go
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except:
    print("You don't seem to have Cython installed. Please get a")
    print("copy from www.cython.org or install it with `pip install Cython`")
    sys.exit(1)


with open('README.md') as f:
    long_description = f.read()

# c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D PRINT_PERSISTENCE_PAIRS
#options = ["-std=c++11", "-Ofast"]
options = ["-std=c++11", "-Ofast"]

if sys.version_info[0] == 2:
    options.append("-fpermissive")


class CustomBuildExtCommand(build_ext):
    """ This extension command lets us not require numpy be installed before running pip install ripser """
    """build_ext command for use when numpy headers are needed."""

    def run(self):
        # Import numpy here, only when headers are needed
        import numpy
        # Add numpy headers to include_dirs
        self.include_dirs.append(numpy.get_include())
        # Call original build_ext command
        build_ext.run(self)


setup(name="ripser",
      version='0.1.7',
      description="A Lean Persistent Homology Library for Python",
      long_description=long_description,
      long_description_content_type="text/markdown",
      author="Chris Tralie, Nathaniel Saul",
      author_email="chris.tralie@gmail.com, nathaniel.saul@wsu.edu",
      url="https://github.com/ctralie/ripser",
      license='LGPL',
      packages=['ripser'],
      ext_modules=cythonize(Extension("pyRipser",
                                      sources=["ripser/pyRipser.pyx"],
                                      define_macros=[("USE_COEFFICIENTS", 1),
                                                     ("NDEBUG", 1), ("ASSEMBLE_REDUCTION_MATRIX", 1)],
                                      extra_compile_args=options,
                                      language="c++"
                                      )),

      install_requires=[
          'Cython',
          'numpy',
          'scipy',
          'matplotlib',
          'scikit-learn'
      ],
      cmdclass={'build_ext': CustomBuildExtCommand},
      )
