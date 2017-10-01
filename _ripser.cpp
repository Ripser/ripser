/*
Programmer: Chris Tralie
Purpose: Python extension/wrapper around ripser
*/

#include <Python.h>
#include <numpy/arrayobject.h>
#include "ripser.h"

/* Docstrings */
static char module_docstring[] =
    "This module provides a wrapper around ripser";
static char ripser_docstring[] =
    "Perform ripser on a distance matrix";

/* Available functions */
static PyObject* ripser_entry(PyObject* self, PyObject* args);

/* Module specification */
static PyMethodDef module_methods[] = {
    {"ripser", ripser_entry, METH_VARARGS, ripser_docstring},
    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef RipserMod =
{
    PyModuleDef_HEAD_INIT,
    "ripser", /* name of module */
    module_docstring, /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};
/* Initialize the module */
PyMODINIT_FUNC PyInit__Ripser(void)
{
    /* Load `numpy` functionality. */
    import_array();
    return PyModule_Create(&RipserMod);
}
#else
/* Initialize the module */
PyMODINIT_FUNC initripser(void)
{
    PyObject *m = Py_InitModule3("ripser", module_methods, module_docstring);
    if (m == NULL)
        return;

    /* Load `numpy` functionality. */
    import_array();
}
#endif

static PyObject *ripser_entry(PyObject *self, PyObject *args)
{
    PyObject *D_obj;
    int modulus, dim_max;
    double threshold;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "Oiid", &D_obj, &modulus, &dim_max, &threshold))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *D_array = PyArray_FROM_OTF(D_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
    if (D_array == NULL) {
        Py_XDECREF(D_array);
        return NULL;
    }

    int N = (int)PyArray_DIM(D_array, 0);

    /* Get pointers to the data as C-types. */
    double* D    = (double*)PyArray_DATA(D_array);

    /* Perform ripser */
    pythondm(D, N, modulus, dim_max, threshold);

    double score = 0.0;

    /* Clean up. */
    Py_DECREF(D_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", score);
    return ret;
}
