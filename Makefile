build: ripser


all: ripser ripser-debug


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG

ripser-debug: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-debug -g


clean:
	rm -f ripser ripser-debug
