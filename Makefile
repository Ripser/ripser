build: ripser


all: ripser ripser-coeff ripser-debug ripser-serial ripser-serial-debug

FLAGS=-Iinclude -pthread
FLAGS+=-MMD -MP

-include *.d

ripser: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o $@ -Ofast -D NDEBUG ${FLAGS}

ripser-serial: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o $@ -Ofast -D NDEBUG ${FLAGS} -DUSE_SERIAL

ripser-serial-debug: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o $@ -g ${FLAGS} -DUSE_SERIAL

ripser-coeff: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o $@ -Ofast -D NDEBUG -D USE_COEFFICIENTS ${FLAGS}

ripser-debug: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o $@ -g ${FLAGS} # -fsanitize=thread -fsanitize=undefined

ripser-tbb: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o $@ -Ofast -D NDEBUG ${FLAGS} -DUSE_TBB -DUSE_TBB_HASHMAP -ltbb

ripser-pstl: ripser.cpp
	c++ -std=c++11 -Wall ripser.cpp -o $@ -Ofast -D NDEBUG ${FLAGS} -DUSE_PARALLEL_STL -DUSE_TBB_HASHMAP -ltbb -std=c++17

clean:
	rm -f ripser ripser-*
