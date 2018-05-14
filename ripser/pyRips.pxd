from libcpp.vector cimport vector

cdef extern from "ripser.cpp":
	vector[float] rips_dm(float* D, int N, int modulus, 
					     int dim_max, float threshold, int do_cocycles)
	
cdef extern from "ripser.cpp":
	vector[float] rips_dm_sparse(int* I, int* J, float* V, int NEdges, 
							   int N, int modulus, int dim_max, 
							   float threshold, int do_cocycles)