"""
    Rips() provides functionality for persistent cohomology calculations, including computing barcodes, cocycles, and visualizations.

"""

from itertools import cycle

import matplotlib.pyplot as plt
import matplotlib as mpl


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
        self.dm_ = None  # Distance matrix
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

        n_points = dm.shape[0]

        if self.thresh == -1:
            thresh = np.max(dm)*2
        else:
            thresh = self.thresh

        [I, J] = np.meshgrid(np.arange(n_points), np.arange(n_points))
        DParam = np.array(dm[I > J], dtype=np.float32)

        res = DRFDM(DParam, self.maxdim, thresh,
                    self.coeff, int(self.do_cocycles))
        pds = []
        for dim in range(self.maxdim + 1):
            # Number of homology classes in this dimension
            n_classes = int(res[0])
            # First extract the persistence diagram
            res = res[1::]
            pd = np.array(res[0:n_classes*2])
            pds.append(np.reshape(pd, (n_classes, 2)))
            res = res[n_classes*2::]
            # Now extract the representative cocycles if they were computed
            if self.do_cocycles and dim > 0:
                self.cocycles_[dim] = []
                for _ in range(n_classes):
                    c_length = int(res[0])
                    res = res[1::]
                    cocycle = res[0:c_length*(dim+2)]
                    cocycle = np.reshape(cocycle, (c_length, dim+2))
                    cocycle = np.array(cocycle, dtype=np.int64)
                    cocycle[:, -1] = np.mod(cocycle[:, -1], self.coeff)
                    self.cocycles_[dim].append(cocycle)
                    res = res[c_length*(dim+2)::]
        return pds

    def plot(self, diagrams=None, plot_only=None, title=None, xy_range=None, labels=None, colormap='default', size=20, ax_color=np.array([0.0, 0.0, 0.0]), colors=None, diagonal=True, lifetime=False, legend=True, show=True):
        """A helper function to plot persistence diagrams. 

        Parameters
        ----------

        diagrams: ndarray (n_pairs, 2) or list of diagrams
            A diagram or list of diagrams as returned from self.fit. If diagram is None, we use self.dgm_ for plotting. If diagram is a list of diagrams, then plot all on the same plot using different colors.

        plot_only: list of numeric
            If specified, an array of only the diagrams that should be plotted.

        title: string, default is None
            If title is defined, add it as title of the plot.

        xy_range: list of numeric [xmin, xmax, ymin, ymax]
            User provided range of axes. This is useful for comparing multiple persistence diagrams.

        labels: string or list of strings
            Legend labels for each diagram. If none are specified, we use H_0, H_1, H_2,... by default.

        colormap: string, default is 'default'
            Any of matplotlib color palettes. Some options are 'default', 'seaborn', 'sequential'. 
            See all availble styles with
            ```
                import matplotlib as mpl
                print(mpl.styles.available)
            ```

        size: numeric, default is 20
            Pixel size of each point plotted.

        ax_color: any valid matplitlib color type. 
            See [https://matplotlib.org/api/colors_api.html](https://matplotlib.org/api/colors_api.html) for complete API.

        diagonal: bool, default is True
            Plot the diagonal x=y line.

        lifetime: bool, default is False. If True, diagonal is turned to False.
            Plot life time of each point instead of birth and death. Essentially, visualize (x, y-x).

        legend: bool, default is True
            If true, show the legend.

        show: bool, default is True
            Call plt.show() after plotting. If you are using self.plot() as part of a subplot, set show=False and call plt.show() only once at the end.

        """

        if labels is None:
            # Provide default labels for diagrams if using self.dgm_
            labels = ["$H_0$", "$H_1$", "$H_2$", "$H_3$",
                      "$H_4$", "$H_5$", "$H_6$", "$H_7$", "$H_8$"]

        if diagrams is None:
            # Allow using transformed diagrams as default
            diagrams = self.dgm_

        if not isinstance(diagrams, list):
            # Must have diagrams as a list for processing downstream
            diagrams = [diagrams]

        if plot_only:
            diagrams = [diagrams[i] for i in plot_only]
            labels = [labels[i] for i in plot_only]

        if not isinstance(labels, list):
            labels = [labels] * len(diagrams)

        if colors is None:
            mpl.style.use(colormap)
            colors = cycle(['C0', 'C1', 'C2', 'C3', 'C4',
                            'C5', 'C6', 'C7', 'C8', 'C9'])

        # Construct copy with proper type of each diagram so we can freely edit them.
        diagrams = [dgm.astype(np.float32, copy=True) for dgm in diagrams]

        # find min and max of all visible diagrams
        concat_dgms = np.concatenate(diagrams).flatten()
        has_inf = np.any(np.isinf(concat_dgms))
        finite_dgms = concat_dgms[np.isfinite(concat_dgms)]

        if not xy_range:
            # define bounds of diagram
            ax_min, ax_max = np.min(finite_dgms), np.max(finite_dgms)
            ax_range = ax_max - ax_min

            # Give plot a nice buffer on all sides.  ax_range=0 when only one point,
            buffer = 1 if ax_range == 0 else ax_range / 5

            ax = ax_min - buffer/2
            bx = ax_max + buffer

            ay, by = ax, bx
        else:
            ax, bx, ay, by = xy_range

        # have inf line slightly below top
        b_inf = bx * 0.95

        xlabel, ylabel = "Birth", "Death"

        if lifetime:
            # Don't plot landscape and diagonal at the same time.
            diagonal = False

            # reset y axis so it doesn't go much below zero
            ay = - buffer/2

            # set custom ylabel
            ylabel = "Lifetime"

            # set diagrams to be (x, y-x)
            for dgm in diagrams:
                dgm[:, 1] -= dgm[:, 0]

            # plot horizon line
            plt.plot([ax, bx], [0, 0], '--', c=ax_color)

        # Plot diagonal
        if diagonal:
            plt.plot([ax, bx], [ax, bx], '--', c=ax_color)

        # Plot inf line
        if has_inf:
            plt.plot([ax, bx], [b_inf, b_inf], c='k', label=r'$\infty$')

            # convert each inf in each diagram with b_inf
            for dgm in diagrams:
                dgm[np.isinf(dgm)] = b_inf

        # Plot each diagram
        for dgm, color, label in zip(diagrams, colors, labels):
            # plot persistence pairs
            plt.scatter(dgm[:, 0], dgm[:, 1], size,
                        color, label=label, edgecolor='none')

            plt.xlabel(xlabel)
            plt.ylabel(ylabel)

        plt.xlim([ax, bx])
        plt.ylim([ay, by])

        if title is not None:
            plt.title(title)

        if legend is True:
            plt.legend(loc='lower right')

        if show is True:
            plt.show()
