build: ripser


all: ripser ripser-coeff ripser-debug

FLAGS=-Iinclude

#FLAGS+=-DUSE_PARALLEL_STL -ltbb -std=c++17

ripser: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser -Ofast -D NDEBUG ${FLAGS}

ripser-coeff: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS ${FLAGS}

ripser-debug: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-debug -g ${FLAGS}


clean:
	rm -f ripser ripser-coeff ripser-debug
