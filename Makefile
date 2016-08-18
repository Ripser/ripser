build: ripser


all: ripser ripser-coeff ripser-reduction


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS

ripser-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-reduction -Ofast -D NDEBUG -D ASSEMBLE_REDUCTION_MATRIX


clean:
	rm -f ripser ripser-coeff ripser-reduction
