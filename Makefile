ifeq ($(platform), Windows)
	EXT := .dll
else
	EXT := .so
endif

build: ripser

all: ripser ripser-coeff ripser-reduction ripser-coeff-reduction ripser-debug libripser$(EXT)

python: src/ripser.cpp
	c++ -std=c++11 src/ripser.cpp -c -o ripser -Ofast -D NDEBUG -D PYTHON_EXTENSION

# ripser: ripser.cpp
# 	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D PRINT_PERSISTENCE_PAIRS

# ripser-coeff: ripser.cpp
# 	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS -D PRINT_PERSISTENCE_PAIRS

# ripser-reduction: ripser.cpp
# 	c++ -std=c++11 ripser.cpp -o ripser-reduction -Ofast -D NDEBUG -D ASSEMBLE_REDUCTION_MATRIX -D PRINT_PERSISTENCE_PAIRS

ripser-coeff-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff-reduction -Ofast -D NDEBUG -D USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX

ripser-debug: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-debug -g

libripser: ripser/ripser.cpp
	c++ -std=c++11 -Ofast -fPIC -shared -L. -D NDEBUG -D USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX -D LIBRIPSER ripser/ripser.cpp -o libripser$(EXT)

clean:
	rm -f ripser ripser-coeff ripser-reduction ripser-coeff-reduction ripser-debug libripser$(EXT)
