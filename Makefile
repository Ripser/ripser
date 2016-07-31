build: ripser


all: ripser ripser-coeff ripser-point-cloud ripser-point-cloud-coeff ripser-dipha ripser-dipha-coeff


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_TRIANGULAR_CSV

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_TRIANGULAR_CSV -D USE_COEFFICIENTS

ripser-point-cloud: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-point-cloud -Ofast -D NDEBUG -D FILE_FORMAT_POINT_CLOUD

ripser-point-cloud-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-point-cloud-coeff -Ofast -D NDEBUG -D FILE_FORMAT_POINT_CLOUD -D USE_COEFFICIENTS

ripser-dipha: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-dipha -Ofast -D NDEBUG -D FILE_FORMAT_DIPHA

ripser-dipha-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-dipha-coeff -Ofast -D NDEBUG -D FILE_FORMAT_DIPHA -D USE_COEFFICIENTS


clean:
	rm ripser ripser-coeff ripser-point-cloud ripser-point-cloud-coeff ripser-dipha ripser-dipha-coeff
