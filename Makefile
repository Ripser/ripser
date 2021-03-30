build: infiltrator


all: infiltrator infiltrator-coeff infiltrator-debug


infiltrator: infiltrator.cpp
	c++ -std=c++11 -Wall infiltrator.cpp -o infiltrator -O3 -D NDEBUG

infiltrator-coeff: infiltrator.cpp
	c++ -std=c++11 -Wall infiltrator.cpp -o infiltrator-coeff -O3 -D NDEBUG -D USE_COEFFICIENTS

infiltrator-debug: infiltrator.cpp
	c++ -std=c++11 -Wall infiltrator.cpp -o infiltrator-debug -g


clean:
	rm -f infiltrator infiltrator-coeff infiltrator-debug
