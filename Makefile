build: ripser-image


all: ripser-image ripser-image-debug


ripser-image: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-image -O3 -D NDEBUG

ripser-image-debug: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-image-debug -g


clean:
	rm -f ripser-image ripser-image-debug
