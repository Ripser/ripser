build: ripser


all: ripser ripser-coeff ripser-reduction ripser-debug


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS

ripser-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-reduction -Ofast -D NDEBUG -D ASSEMBLE_REDUCTION_MATRIX

ripser-debug: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-debug -g


clean:
	rm -f ripser ripser-coeff ripser-reduction ripser-debug


live: index.html common.js ripser-web.js ripser-worker.js 
	source build_emscripten.sh && NACL_SDK_ROOT=~/Source/nacl_sdk/pepper_49/ make -f Makefile_pnacl && cp -r index.html common.js ripser-web.js ripser-worker.js pnacl emscripten ../
