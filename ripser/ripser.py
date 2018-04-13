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

    Example usage:

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

    def transform(self, X, distance_matrix=False, metric='euclidean'):

        if not distance_matrix:
            X = pairwise_distances(X, metric=metric)

        dgm = self.compute_rips(X)
        self._dgm = dgm
    
    def fit_transform(self, X, distance_matrix=False, metric='euclidean'):
        """ Run transform and return results

        """
        self.transform(X, distance_matrix, metric)
        return self._dgm

    def compute_rips(self, dm):
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

    def plot(self, diagram=None, diagonal=True, sz=20, label='dgm', axcolor=np.array([0.0, 0.0, 0.0]), marker=None, show=True):
        """ Plot each diagram on the same plot.
        """
        if diagram is None:
            diagram = self._dgm
        
        if type(diagram) is not list:
            diagram = [diagram]
    
        colors = ['r','g', 'b'] # TODO: convert this to a cylic generator so we can zip as many as required.
        for dgm, color in zip(diagram, colors):
            
            if dgm.size is not 0:
                # build diagonal line

                if diagonal:
                    axMin, axMax = np.min(dgm), np.max(dgm)
                    axRange = axMax - axMin

                    a = max(axMin - axRange / 5, 0)
                    b = axMax + axRange / 5

                    # plot diagonal line
                    plt.plot([a, b], [a, b], '--', c=axcolor, label='none')
                
                # plot points
                plt.scatter(dgm[:, 0], dgm[:, 1], sz, 
                                color, label=label, edgecolor='none')
                # add labels
                plt.xlabel('Birth')
                plt.ylabel('Death')
                # plt.axis('equal')

        if show:
            plt.show()



