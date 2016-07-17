all:
	c++ -std=c++11 ripser.cpp -o ripser -Ofast -D NDEBUG -D FILE_FORMAT_LOWER_TRIANGULAR_CSV -D PRINT_PERSISTENCE_PAIRS
