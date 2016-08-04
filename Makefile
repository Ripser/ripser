build: ripser


all: ripser ripser-coeff ripser-dist ripser-dist-coeff ripser-upper-dist ripser-upper-dist-coeff ripser-point-cloud ripser-point-cloud-coeff ripser-dipha ripser-dipha-coeff


ripser: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_DISTANCE_MATRIX

ripser-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_DISTANCE_MATRIX -D USE_COEFFICIENTS

ripser-dist: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-dist -Ofast -D NDEBUG -D FILE_FORMAT_DISTANCE_MATRIX

ripser-dist-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-dist-coeff -Ofast -D NDEBUG -D FILE_FORMAT_DISTANCE_MATRIX -D USE_COEFFICIENTS

ripser-upper-dist: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-upper-dist -Ofast -D NDEBUG -D FILE_FORMAT_UPPER_DISTANCE_MATRIX

ripser-upper-dist-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-upper-dist-coeff -Ofast -D NDEBUG -D FILE_FORMAT_UPPER_DISTANCE_MATRIX -D USE_COEFFICIENTS

ripser-point-cloud: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-point-cloud -Ofast -D NDEBUG -D FILE_FORMAT_POINT_CLOUD

ripser-point-cloud-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-point-cloud-coeff -Ofast -D NDEBUG -D FILE_FORMAT_POINT_CLOUD -D USE_COEFFICIENTS

ripser-dipha: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-dipha -Ofast -D NDEBUG -D FILE_FORMAT_DIPHA

ripser-dipha-coeff: ripser.cpp
	c++ -std=c++11 ripser.cpp -o ripser-dipha-coeff -Ofast -D NDEBUG -D FILE_FORMAT_DIPHA -D USE_COEFFICIENTS


clean:
	rm ripser ripser-coeff ripser-dist ripser-dist-coeff ripser-upper-dist ripser-upper-dist-coeff ripser-point-cloud ripser-point-cloud-coeff ripser-dipha ripser-dipha-coeff
