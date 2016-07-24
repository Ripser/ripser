build: ripser


all: ripser ripser_coeff ripser_dipha ripser_dipha_coeff


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_TRIANGULAR_CSV

ripser_coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser_coeff -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_TRIANGULAR_CSV -D USE_COEFFICIENTS

ripser_dipha: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser_dipha -Ofast -D NDEBUG -D FILE_FORMAT_DIPHA

ripser_dipha_coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser_dipha_coeff -Ofast -D NDEBUG -D FILE_FORMAT_DIPHA -D USE_COEFFICIENTS


clean:
	rm ripser ripser_coeff ripser_dipha ripser_dipha_coeff
