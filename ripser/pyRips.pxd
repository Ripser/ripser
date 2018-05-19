from libcpp.vector cimport vector

cdef extern from "ripser.cpp":
	ctypedef struct ripserResults:
		vector[vector[float]] births_and_deaths_by_dim
		vector[vector[vector[int]]] cocycles_by_dim
		int num_edges

cdef extern from "ripser.cpp":
	ripserResults rips_dm(float* D, int N, int modulus, 
					     int dim_max, float threshold, int do_cocycles)
	
cdef extern from "ripser.cpp":
	ripserResults rips_dm_sparse(int* I, int* J, float* V, int NEdges, 
							   int N, int modulus, int dim_max, 
							   float threshold, int do_cocycles)