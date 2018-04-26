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
    do_cocycles: bool
        Indicator of whether to compute cocycles, if so, we compute and store
        cocycles in the cocycles_ dictionary Rips member variable
    
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

    def __init__(self, maxdim=1, thresh=-1, coeff=2, do_cocycles=False, verbose=True):
        self.maxdim = maxdim
        self.thresh = thresh
        self.coeff = coeff
        self.do_cocycles = do_cocycles
        self.verbose = verbose

        self.dgm_ = None
        self.cocycles_ = {}
        self.dm_ = None #Distance matrix
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
        X = np.array(X, dtype=np.float32)
        self.dm_ = X
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
        """
        
        npts = dm.shape[0]
        if self.thresh == -1:
            thresh = np.max(dm)*2
        else:
            thresh = self.thresh
        [I, J] = np.meshgrid(np.arange(npts), np.arange(npts))
        DParam = np.array(dm[I > J], dtype=np.float32)

        res = DRFDM(DParam, self.maxdim, thresh, self.coeff, int(self.do_cocycles))
        pds = []
        for dim in range(self.maxdim + 1):
            nclasses = int(res[0]) #Number of homology classes in this dimension
            #First extract the persistence diagram
            res = res[1::]
            pd = np.array(res[0:nclasses*2])
            pds.append(np.reshape(pd, (nclasses, 2)))
            res = res[nclasses*2::]
            #Now extract the representative cocycles if they were computed
            if self.do_cocycles and dim > 0:
                self.cocycles_[dim] = []
                for n in range(nclasses):
                    clen = int(res[0])
                    res = res[1::]
                    cocycle = res[0:clen*(dim+2)]
                    cocycle = np.reshape(cocycle, (clen, dim+2))
                    cocycle = np.array(cocycle, dtype=np.int64)
                    cocycle[:, -1] = np.mod(cocycle[:, -1], self.coeff)
                    self.cocycles_[dim].append(cocycle)
                    res = res[clen*(dim+2)::]
        return pds

    def plot(self, diagrams=None, plotonly=None, diagonal=True, sz=20, labels=None, axcolor=np.array([0.0, 0.0, 0.0]), colors=None, marker=None, title=None, legend=True, show=True):
        """A helper function to plot persistence diagrams.

        Parameters
        ----------
        diagrams: ndarray (n_pairs, 2) or list of diagrams
            A diagram or list of diagrams as returned from self.fit. If diagram is None, we use self.dgm_ for plotting. If diagram is a list of diagrams, then plot all on the same plot using different colors.
        
        plotonly: If specified, an array of only the diagrams that should be plotted

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


        # find min and max of all visible diagrams, plot diagonal only once
        concatdgms = np.array([])
        if not plotonly:
            plotonly = range(len(diagrams))
        concatdgms = np.concatenate([diagrams[i] for i in plotonly]).flatten()
        numinf = np.sum(np.isinf(concatdgms))
        concatdgms[np.isinf(concatdgms)] = np.min(concatdgms)
        axMin, axMax = np.min(concatdgms), np.max(concatdgms)
        axRange = axMax - axMin
        a = max(axMin - axRange / 5, 0)
        b = axMax + axRange / 5
        if a == b:
            a = 0
            b = 1
        fuzz = 0.05*(b-a)
        a -= fuzz

        if diagonal:
            plt.plot([a, b], [a, b], '--', c=axcolor)
        
        if numinf > 0:
            plt.plot([a-fuzz, b+fuzz], [b*1.05]*2, c='k')
            plt.text(a-fuzz, b+fuzz, '$\infty$', size=14)

        # Plot each diagram
        for i, (dgm, color, label) in enumerate(zip(diagrams, colors, labels)):
            if dgm.size is not 0 and i in plotonly:
                # plot persistence pairs
                finitedgm = dgm[np.isfinite(dgm[:, 1]), :]
                plt.scatter(finitedgm[:, 0], finitedgm[:, 1], sz, 
                                color, label=label, edgecolor='none')
                infdgm = np.array(dgm[np.isinf(dgm[:, 1]), :])
                infdgm[:, 1] = b*1.1
                plt.scatter(infdgm[:, 0], infdgm[:, 1], sz, 
                                color, edgecolor='none')
                plt.xlabel('Birth')
                plt.ylabel('Death')

        
        plt.xlim([a-fuzz, b+fuzz])
        plt.ylim([a-fuzz, b*1.1+fuzz])

        if title is not None:
            plt.title(title)

        if legend is True:
            plt.legend(loc='lower right')

        if show is True:
            plt.show()



