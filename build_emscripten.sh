source ~/Source/emsdk_portable/emsdk_env.sh
mkdir -p emscripten
emcc -s ALLOW_MEMORY_GROWTH=1 -Wall --bind --memory-init-file 0 -std=c++11 -o emscripten/ripser.js -O3 -D NDEBUG ripser.cpp

