from libcpp.vector cimport vector

cdef extern from "ripser.cpp":
	vector[float] pythondm(float* D, int N, int modulus, 
					int dim_max, float threshold)
	
