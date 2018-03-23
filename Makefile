build: ripser-representatives


ripser-representatives: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-representatives -Ofast -g -D NDEBUG -D _USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS

ripser-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-reduction -Ofast -D NDEBUG -D ASSEMBLE_REDUCTION_MATRIX

ripser-coeff-reduction: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff-reduction -Ofast -D NDEBUG -D USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX

ripser-debug: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-debug -g


clean:
	rm -f ripser-representatives
