# Ripser

Copyright © 2015–2016 [Ulrich Bauer].


### Description

Ripser is a lean C++ code for the computation of Vietoris–Rips persistence barcodes. It can do just this one thing, but does it extremely well.

The main features of Ripser:

  - time- and memory-efficient
  - less than 1000 lines of code in a single C++ file
  - support for coefficients in prime finite fields
  - no external dependencies (optional support for Google's [sparsehash])

Currently, Ripser outperforms other codes ([Dionysus], [DIPHA], [GUDHI], [Perseus], [PHAT]) by a factor of more than 40 in computation time and a factor of more than 15 in memory efficiency. (Note that [PHAT] does not contain code for generating Vietoris–Rips filtrations).

Input formats currently supported by Ripser:

  - comma-separated values lower triangular distance matrix (preferred)
  - comma-separated values upper triangular distance matrix (MATLAB output from the function `pdist`)
  - [DIPHA] distance matrix data

Ripser's efficiency is based on a few important concepts and principles:
  
  - Compute persistent *co*homology
  - Don't compute information that is never needed
    (for the experts: employ the *clearing* optimization, aka *persistence with a twist*)
  - Don't store information that can be readily recomputed
  - Take obvious shortcuts (*apparent persistence pairs*)


### Version
[Latest release][latest-release]: 1.0 (July 2016)


### Building

Ripser requires a C++11 compiler. Here is how to obtain, build, and run Ripser:

```sh
git clone https://github.com/Ripser/ripser.git
cd ripser
make
./ripser examples/sphere_3_192.lower_distance_matrix
```


### Options

Ripser supports several compile-time options. They are switched on by defining the C preprocessor macros listed below, either using `#define` in the code or by passing an argument to the compiler. The following options are supported:

  - `ASSEMBLE_REDUCTION_MATRIX`: store the reduction matrix; may speed up computation but will also increase memory usage
  - `USE_COEFFICIENTS`: enable support for coefficients in a prime field
  - `INDICATE_PROGRESS`: indicate the current progress in the console
  - `PRINT_PERSISTENCE_PAIRS`: output the computed persistence pairs (enabled by default)
  - `USE_GOOGLE_HASHMAP`: enable support for Google's [sparsehash] data structure; may further reducue memory footprint

Furthermore, one of the following options needs to be chosen to specify the format for the input files:

  - `FILE_FORMAT_LOWER_TRIANGULAR_CSV`: lower triangular distance matrix; a comma (or whitespace, or other non-numerical character) separated list of the distance matrix entries below the diagonal, sorted lexicographically by row index, then column index
  - `FILE_FORMAT_UPPER_TRIANGULAR_CSV`: upper triangular distance matrix; similar to the previous, but for the entries above the diagonal; suitable for output from the MATLAB function `pdist`, saved in a CSV file
  - `FILE_FORMAT_DIPHA`: DIPHA distance matrix as described on the [DIPHA] website

For example, to build with support for coefficients:

```sh
$ c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_TRIANGULAR_CSV -D USE_COEFFICIENTS
```

The following options are supported at the command line:

  - `--dim k`: compute persistent homology up to dimension *k*
  - `--threshold t`: compute Rips complexes up to diameter *t*
  - `--modulus p`: compute homology with coefficients in the prime field Z/*p*Z (only available when build with the option `USE_COEFFICIENTS`)


### Planned features

The following features are currently planned for future versions:

 - support for point clouds
 - computation of representative cycles for persistent homology (currenly only *co*cycles are computed)
 - support for sparse distance matrices


### License

[LGPL] 3.0


[Ulrich Bauer]: <http://ulrich-bauer.org>
[latest-release]: <https://github.com/Ripser/ripser/releases/latest>
[Dionysus]: <http://www.mrzv.org/software/dionysus/>
[DIPHA]: <http://git.io/dipha>
[PHAT]: <http://git.io/dipha>
[Perseus]: <http://www.sas.upenn.edu/~vnanda/perseus/>
[GUDHI]: <http://gudhi.gforge.inria.fr>
[sparsehash]: <https://github.com/sparsehash/sparsehash>
[LGPL]: <https://www.gnu.org/licenses/lgpl>