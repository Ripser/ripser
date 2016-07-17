build:
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_TRIANGULAR_CSV

coefficients:
	c++ -std=c++11 ripser.cpp -o ripser_coeff -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_TRIANGULAR_CSV -D USE_COEFFICIENTS
