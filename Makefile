build: ripser


all: ripser ripser-coeff ripser-reduction ripser-debug


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D PRINT_PERSISTENCE_PAIRS

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS -D PRINT_PERSISTENCE_PAIRS

ripser-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-reduction -Ofast -D NDEBUG -D ASSEMBLE_REDUCTION_MATRIX -D PRINT_PERSISTENCE_PAIRS

ripser-debug: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-debug -g


clean:
	rm -f ripser ripser-coeff ripser-reduction ripser-debug
