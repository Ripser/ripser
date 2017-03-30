build: ripser-representatives


ripser-representatives: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-representatives -Ofast -D NDEBUG -D USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX -D USE_GOOGLE_HASHMAP

clean:
	rm -f ripser-representatives
