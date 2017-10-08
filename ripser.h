#include <Python.h>
#include <numpy/arrayobject.h>
std::vector<float> pythondm(float* D, int N, int modulus, int dim_max, float threshold);
