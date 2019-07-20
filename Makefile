build: ripser


all: ripser ripser-coeff ripser-debug


ripser: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser -Ofast -D NDEBUG

ripser-coeff: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS

ripser-debug: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-debug -g


clean:
	rm -rf ripser ripser-coeff ripser-debug emscripten pnacl && NACL_SDK_ROOT=~/Source/nacl_sdk/pepper_49/ make -f Makefile_pnacl clean


live: pnacl emscripten
	cp -r pnacl emscripten ../

pnacl: ripser.cpp
	NACL_SDK_ROOT=~/Source/nacl_sdk/pepper_49/ make -f Makefile_pnacl

emscripten: ripser.cpp
	. build_emscripten.sh
