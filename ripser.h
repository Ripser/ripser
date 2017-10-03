#include <Python.h>
#include <numpy/arrayobject.h>
PyArrayObject* pythondm(float* D, int N, int modulus, int dim_max, float threshold);
