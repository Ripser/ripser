#include <Python.h>
#include <numpy/arrayobject.h>
PyArrayObject* pythondm(double* D, int N, int modulus, int dim_max, double threshold);
