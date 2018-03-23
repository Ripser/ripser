build: ripser-representatives


ripser-representatives: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-representatives -Ofast -g -D NDEBUG -D USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX

clean:
	rm -f ripser-representatives
