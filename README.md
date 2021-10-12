# Ripser for image persistence

Copyright © 2015–2021 [Ulrich Bauer], [Maximilian Schmahl].


### Description

Ripser is a lean C++ code for the computation of Vietoris–Rips persistence barcodes. It can do just this one thing, but does it extremely well.

This branch computes the image persistence barcode induced by a nonexpanding map of metric spaces with the same underlying set, specified by two distance matrices of the same size.

The input is given either in a file whose name is passed as an argument, or through stdin. The following options are supported at the command line:

  - `--format`: use the specified file format for the input.  The following formats are supported:
    - `lower-distance` (default if no format is specified): lower triangular distance matrix; a comma (or whitespace, or other non-numerical character) separated list of the distance matrix entries below the diagonal, sorted lexicographically by row index, then column index
    - `distance`: full distance matrix; similar to the above, but for all entries of the distance matrix
    - `point-cloud`: point cloud; a comma (or whitespace, or other non-numerical character)  separated list of coordinates of the points in some Euclidean space, one point per line
    - `binary`: lower distance matrix in binary file format; a sequence of the distance matrix entries below the diagonal in 64 bit double format (IEEE 754, little endian).
  - `--subfiltration <f>`: use f as second filtration for image persistence
  - `--dim k`: compute persistent homology up to dimension *k*
  - `--threshold t`: compute Rips complexes up to diameter *t*





### License

Ripser for image persistence is licensed under the [LGPL] 3.0. Please contact the authors if you want to use Ripser in your software under a different license. 

[Ulrich Bauer]: <http://ulrich-bauer.org>
[Maximilian Schmahl]: <mailto:maximilian.schmahl@web.de>
[live.ripser.org]: <http://live.ripser.org>
[PNaCl]: <https://www.chromium.org/nativeclient/pnacl/>
[Emscripten]: <http://emscripten.org>
[latest-release]: <https://github.com/Ripser/ripser/releases/latest>
[Dionysus]: <http://www.mrzv.org/software/dionysus/>
[DIPHA]: <http://git.io/dipha>
[PHAT]: <http://git.io/dipha>
[Perseus]: <http://www.sas.upenn.edu/~vnanda/perseus/>
[GUDHI]: <http://gudhi.gforge.inria.fr>
[sparsehash]: <https://github.com/sparsehash/sparsehash>
[LGPL]: <https://www.gnu.org/licenses/lgpl>
