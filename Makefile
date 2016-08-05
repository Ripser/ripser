build: ripser


all: ripser ripser-coeff


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS


clean:
	rm ripser ripser-coeff