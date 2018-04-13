# distutils: language = c++

cimport numpy as np
cimport ripser.pyRips as pyRips
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
def doRipsFiltrationDM(np.ndarray[float,ndim=1,mode="c"] DParam not None, int maxHomDim, int thresh=-1, int coeff=2):

	cdef int N = DParam.shape[0]

	res = pyRips.pythondm(&DParam[0], N, coeff, maxHomDim, thresh)

	return res
