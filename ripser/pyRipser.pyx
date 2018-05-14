# distutils: language = c++

cimport numpy as np
cimport ripser.pyRips as pyRips
import cython

@cython.boundscheck(False)
@cython.wraparound(False)
def doRipsFiltrationDM(np.ndarray[float,ndim=1,mode="c"] DParam not None, int maxHomDim, float thresh=-1, int coeff=2, int do_cocycles=0):

	cdef int N = DParam.shape[0]

	res = pyRips.rips_dm(&DParam[0], N, coeff, maxHomDim, thresh, do_cocycles)

	return res

@cython.boundscheck(False)
@cython.wraparound(False)
def doRipsFiltrationDMSparse(np.ndarray[int,ndim=1,mode="c"] I not None, np.ndarray[int,ndim=1,mode="c"] J not None, np.ndarray[float,ndim=1,mode="c"] V not None, int N, int maxHomDim, float thresh=-1, int coeff=2, int do_cocycles=0):

	cdef int NEdges = I.size

	res = pyRips.rips_dm_sparse(&I[0], &J[0], &V[0], NEdges, N, coeff, maxHomDim, thresh, do_cocycles)

	return res