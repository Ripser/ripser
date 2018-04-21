import subprocess
import os

import time

import matplotlib.pyplot as plt

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.metrics.pairwise import pairwise_distances

from pyRipser import doRipsFiltrationDM as DRFDM


class Rips(BaseEstimator):
    """Wrapper around Uli Bauer's Ripser code 


    Parameters
    ----------
    maxdim : int, optional, default 1
        Maximum homology dimension computed. Will compute all dimensions lower than
        and equal to this value. For 1, H_0 and H_1 will be computed.
    thresh : float, default -1
        Maximum distances considered when constructing filtration. If -1, compute 
        the entire filtration.
    coeff : int prime, default 2
        Compute homology with coefficients in the prime field Z/pZ for p=coeff.
    
    
    Attributes
    ----------
    _dgm : list of ndarray, each shape (n_pairs, 2)
        After `transform`, _dgm contains computed persistence diagrams in
        each dimension

    Examples
    --------

    ```
    from ripser import Rips
    from sklearn import datasets

    data = datasets.make_circles(n_samples=110)[0]
    rips = Rips()
    rips.transform(data)
    rips.plot()
    ```

    """ 



    def __init__(self, maxdim=1, thresh=-1, coeff=2):
        self.maxdim = maxdim
        self.thresh = thresh
        self.coeff = coeff
        self._dgm = None

    def transform(self, X, distance_matrix=False, metric='euclidean'):
        """Compute persistence diagrams for X data array.

        Parameters
        ----------
        X: ndarray (n_samples, n_features)
            A numpy array of either data or distance matrix.
        distance_matrix: bool
            Indicator that X is a distance matrix, if not we compute a 
            distance matrix from X using the chosen metric.
        metric: string or callable
            The metric to use when calculating distance between instances in a feature array. If metric is a string, it must be one of the options specified in PAIRED_DISTANCES, including “euclidean”, “manhattan”, or “cosine”. Alternatively, if metric is a callable function, it is called on each pair of instances (rows) and the resulting value recorded. The callable should take two arrays from X as input and return a value indicating the distance between them.

        """


        # Default is to input point cloud data
        if not distance_matrix:
            X = pairwise_distances(X, metric=metric)

        dgm = self._compute_rips(X)
        self._dgm = dgm

        return dgm
    
    def fit_transform(self, X, distance_matrix=False, metric='euclidean'):
        """ Run transform and return results

        """
        self.transform(X, distance_matrix, metric)
        return self._dgm

    def _compute_rips(self, dm):
        """ Compute the persistence diagram
       
        :param D: An NxN pairwise distance matrix
        :param maxHomDim: The dimension up to which to compute persistent homology
        :param thresh: Threshold up to which to add edges.  If not specified, add all
            edges up to the full clique
        :param coeff: A prime to use as the field coefficients for the PH computation

        :return: PDs (array of all persistence diagrams from 0D up to maxHomDim).
            Each persistence diagram is a numpy array
        """
        
        N = dm.shape[0]
        if self.thresh == -1:
            thresh = np.max(dm)*2
        else:
            thresh = self.thresh
        [I, J] = np.meshgrid(np.arange(N), np.arange(N))
        DParam = np.array(dm[I > J], dtype=np.float32)

        res = DRFDM(DParam, self.maxdim, thresh, self.coeff)
        PDs = []
        istart = 0
        for dim in range(self.maxdim + 1):
            N = int(res[-1 - (self.maxdim - dim)])
            if dim > 0:
                N -= int(res[-2 - (self.maxdim - dim)])
            I = np.array(res[istart * 2:(istart + N) * 2])
            PDs.append(np.reshape(I, (N, 2)))
            istart += N
        return PDs

    def plot(self, diagram=None, diagonal=True, sz=20, labels='dgm', axcolor=np.array([0.0, 0.0, 0.0]), marker=None, show=True):
        """ Plot each diagram on the same plot.
        """
        if diagram is None:
            diagram = self._dgm
        
        if type(diagram) is not list:
            diagram = [diagram]

        
        if type(labels) is not list:
            labels = [labels] * len(diagram)
    
        colors = ['r','g', 'b'] # TODO: convert this to a cylic generator so we can zip as many as required.
        for dgm, color, label in zip(diagram, colors, labels):
            
            if dgm.size is not 0:
                # build diagonal line

                if diagonal:
                    axMin, axMax = np.min(dgm), np.max(dgm)
                    axRange = axMax - axMin

                    a = max(axMin - axRange / 5, 0)
                    b = axMax + axRange / 5

                    # plot diagonal line
                    plt.plot([a, b], [a, b], '--', c=axcolor)
                
                # plot points
                plt.scatter(dgm[:, 0], dgm[:, 1], sz, 
                                color, label=label, edgecolor='none')
                # add labels
                plt.xlabel('Birth')
                plt.ylabel('Death')
                # plt.axis('equal')

        if show:
            plt.show()



