build: ripser-representatives


ripser-representatives: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o ripser-representatives -Ofast -D NDEBUG

clean:
	rm -f ripser-representatives
