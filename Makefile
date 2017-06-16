build: ripser


all: ripser ripser-coeff ripser-reduction ripser-coeff-reduction ripser-debug


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D PRINT_PERSISTENCE_PAIRS -D INDICATE_PROGRESS 

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS

ripser-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-reduction -Ofast -D NDEBUG -D ASSEMBLE_REDUCTION_MATRIX

ripser-coeff-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff-reduction -Ofast -D NDEBUG -D USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX

ripser-debug: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-debug -g


clean:
	yes| rm ripser ripser-coeff ripser-reduction ripser-coeff-reduction ripser-debug emscripten pnacl 
    #&& NACL_SDK_ROOT=~/Source/nacl_sdk/pepper_49/ make -f Makefile_pnacl clean


live: pnacl emscripten
	cp -r index.html common.js ripser-web.js ripser-worker.js pnacl emscripten ../

pnacl: index.html common.js ripser-web.js ripser-worker.js
	NACL_SDK_ROOT=~/Source/nacl_sdk/pepper_49/ make -f Makefile_pnacl

emscripten: index.html common.js ripser-web.js ripser-worker.js
	. build_emscripten.sh
