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
    dgm_ : list of ndarray, each shape (n_pairs, 2)
        After `transform`, dgm_ contains computed persistence diagrams in
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

    def __init__(self, maxdim=1, thresh=-1, coeff=2, verbose=True):
        self.maxdim = maxdim
        self.thresh = thresh
        self.coeff = coeff
        self.verbose = verbose

        self.dgm_ = None
        self.distance_matrix_ = None # indicator
        self.metric_ = None

        if self.verbose:
            print("Rips(maxdim={}, thres={}, coef={}, verbose={})".format(
                maxdim, thresh, coeff, verbose))

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

        if not distance_matrix:
            X = pairwise_distances(X, metric=metric)

        dgm = self._compute_rips(X)
        self.dgm_ = dgm
    
    def fit_transform(self, X, distance_matrix=False, metric='euclidean'):
        """Compute persistence diagrams for X data array and return the diagrams.

        Parameters
        ----------
        X: ndarray (n_samples, n_features)
            A numpy array of either data or distance matrix.
        
        distance_matrix: bool
            Indicator that X is a distance matrix, if not we compute a 
            distance matrix from X using the chosen metric.
        
        metric: string or callable
            The metric to use when calculating distance between instances in a feature array. If metric is a string, it must be one of the options specified in PAIRED_DISTANCES, including “euclidean”, “manhattan”, or “cosine”. Alternatively, if metric is a callable function, it is called on each pair of instances (rows) and the resulting value recorded. The callable should take two arrays from X as input and return a value indicating the distance between them.
        
        Return
        ------
        dgms: list (size maxdim) of ndarray (n_pairs, 2)
            A list of persistence diagrams, one for each dimension less than maxdim. Each diagram is an ndarray of size (n_pairs, 2) with the first column representing the birth time and the second column representing the death time of each pair.

        """
        self.transform(X, distance_matrix, metric)
        return self.dgm_


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

    def plot(self, diagrams=None, diagonal=True, sz=20, labels=None, axcolor=np.array([0.0, 0.0, 0.0]), colors=None, marker=None, title=None, legend=True, show=True):
        """A helper function to plot persistence diagrams.

        Parameters
        ----------
        diagrams: ndarray (n_pairs, 2) or list of diagrams
            A diagram or list of diagrams as returned from self.fit. If diagram is None, we use self.dgm_ for plotting. If diagram is a list of diagrams, then plot all on the same plot using different colors.
        
        diagonal: bool, default is True
            Plot the diagonal line
        
        labels: string or list of strings
            Legend labels for each diagram. If none are specified, assume the first diagram is H_0 and we move up from there.
        
        title: string, default is None
            If title is defined, add it as title of the plot.

        legend: bool, default is True
            If true, show the legend.

        show: bool, default is True
            Call plt.show() after plotting. If you are using self.plot() as part of a subplot, set show=False and call plt.show() only once at the end.
        
        """

        ##  Note:  This method is a bit overloaded to accomodate a 
        #          single diagram or a list of diagrams. Please keep this 
        #          in mind when making changes. 
        #          Refactors are welcome.

        if labels is None:
            # Provide default labels for diagrams if using self.dgm_
            labels = ["H0", "H1", "H2", "H3", "H4", "H5", "H6", "H8"]
        if diagrams is None:
            # Allow using transformed diagrams as default
            diagrams = self.dgm_
        if type(diagrams) is not list:
            # Must have diagrams as a list
            diagrams = [diagrams]
        if type(labels) is not list:
            labels = [labels] * len(diagrams)
        if colors is None:
            # TODO: convert this to a cylic generator so we can zip as many as required.
            colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'] 


        # find min and max of all diagrams, plot diagonal only once
        if diagonal:
            axMin, axMax = np.min(np.concatenate(diagrams)), np.max(np.concatenate(diagrams))
            axRange = axMax - axMin

            a = max(axMin - axRange / 5, 0)
            b = axMax + axRange / 5
            plt.plot([a, b], [a, b], '--', c=axcolor)

        # Plot each diagram
        for dgm, color, label in zip(diagrams, colors, labels):
            if dgm.size is not 0:
                # plot persistence pairs
                plt.scatter(dgm[:, 0], dgm[:, 1], sz, 
                                color, label=label, edgecolor='none')
                
                plt.xlabel('Birth')
                plt.ylabel('Death')

        if title is not None:
            plt.title(title)

        if legend is True:
            plt.legend(loc='lower right')

        if show is True:
            plt.show()



