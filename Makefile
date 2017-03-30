build: ripser


all: ripser ripser-coeff ripser-reduction ripser-coeff-reduction ripser-debug


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D USE_GOOGLE_HASHMAP

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_GOOGLE_HASHMAP -D USE_COEFFICIENTS

ripser-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-reduction -Ofast -D NDEBUG -D USE_GOOGLE_HASHMAP -D ASSEMBLE_REDUCTION_MATRIX

ripser-coeff-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff-reduction -Ofast -D NDEBUG -D USE_GOOGLE_HASHMAP -D USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX

ripser-debug: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-debug -g -D USE_GOOGLE_HASHMAP


clean:
	rm -f ripser ripser-coeff ripser-reduction ripser-coeff-reduction ripser-debug
