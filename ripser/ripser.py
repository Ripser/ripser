"""
    ripser provides functionality for persistent cohomology calculations, 
    including computing barcodes, cocycles, and visualizations.

"""

from itertools import cycle
import warnings

import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import sparse

import numpy as np
from sklearn.base import TransformerMixin
from sklearn.metrics.pairwise import pairwise_distances

from pyRipser import doRipsFiltrationDM as DRFDM
from pyRipser import doRipsFiltrationDMSparse as DRFDMSparse


def ripser(X, maxdim=1, thresh=np.inf, coeff=2, distance_matrix=False,
           do_cocycles=False, metric='euclidean'):
    """ Compute persistence diagrams for X data array. If X is not a 
        distance matrix, it will be converted to a distance matrix using 
        the chosen metric.

    Parameters
    ----------
    X: ndarray (n_samples, n_features)
        A numpy array of either data or distance matrix.
        Can also be a sparse distance matrix of type scipy.sparse

    maxdim : int, optional, default 1
        Maximum homology dimension computed. Will compute all dimensions 
        lower than and equal to this value. 
        For 1, H_0 and H_1 will be computed.

    thresh : float, default infinity
        Maximum distances considered when constructing filtration. 
        If infinity, compute the entire filtration.

    coeff : int prime, default 2
        Compute homology with coefficients in the prime field Z/pZ for p=coeff.

    distance_matrix: bool
        Indicator that X is a distance matrix, if not we compute a 
        distance matrix from X using the chosen metric.

    do_cocycles: bool
        Indicator of whether to compute cocycles, if so, we compute and store
        cocycles in the cocycles_ dictionary Rips member variable

    metric: string or callable
        The metric to use when calculating distance between instances in a 
        feature array. If metric is a string, it must be one of the options 
        specified in PAIRED_DISTANCES, including "euclidean", "manhattan", 
        or "cosine". Alternatively, if metric is a callable function, it is 
        called on each pair of instances (rows) and the resulting value 
        recorded. The callable should take two arrays from X as input and 
        return a value indicating the distance between them.

    Return
    ------
    A dictionary holding all of the results of the computation

    {'dgms': list (size maxdim) of ndarray (n_pairs, 2)
        A list of persistence diagrams, one for each dimension less 
        than maxdim. Each diagram is an ndarray of size (n_pairs, 2) 
        with the first column representing the birth time and the 
        second column representing the death time of each pair.
     'cocycles': list (size maxdim)
        A list of representative cocycles in each dimension.  The list 
        in each dimension is parallel to the diagram in that dimension.
     'num_edges': int
        The number of edges added during the computation
     'dm' : ndarray (n_samples, n_samples)
        The distance matrix used in the computation
    }

    Examples
    --------

    ```
    from ripser import ripser, plot_dgms
    from sklearn import datasets

    data = datasets.make_circles(n_samples=110)[0]
    dgms = ripser(data)['dgms']
    plot_dgms(dgms)
    ```

    """

    if not distance_matrix:
        if X.shape[0] == X.shape[1]:
            warnings.warn(
                "The input matrix is square, but the distance_matrix " +
                "flag is off.  Did you mean to indicate that " +
                "this was a distance matrix?")
        elif X.shape[0] < X.shape[1]:
            warnings.warn(
                "The input point cloud has more columns than rows; " +
                "did you mean to transpose?")
        X = pairwise_distances(X, metric=metric)

    if not (X.shape[0] == X.shape[1]):
        raise Exception('Distance matrix is not square')
    dm = X
    n_points = dm.shape[0]

    if sparse.issparse(dm):
        coo = sparse.coo_matrix.astype(dm.tocoo(), dtype=np.float32)
        res = DRFDMSparse(coo.row, coo.col, coo.data, n_points,
                          maxdim, thresh, coeff, int(do_cocycles))
    else:
        I, J = np.meshgrid(np.arange(n_points), np.arange(n_points))
        DParam = np.array(dm[I > J], dtype=np.float32)
        res = DRFDM(DParam, maxdim, thresh, coeff, int(do_cocycles))

    # Unwrap persistence diagrams
    dgms = res['births_and_deaths_by_dim']
    for dim in range(len(dgms)):
        N = int(len(dgms[dim])/2)
        dgms[dim] = np.reshape(np.array(dgms[dim]), [N, 2])

    # Unwrap cocycles
    cocycles = []
    for dim in range(len(res['cocycles_by_dim'])):
        cocycles.append([])
        for j in range(len(res['cocycles_by_dim'][dim])):
            ccl = res['cocycles_by_dim'][dim][j]
            n = int(len(ccl)/(dim+2))
            ccl = np.reshape(np.array(ccl, dtype=np.int64), [n, dim+2])
            ccl[:, -1] = np.mod(ccl[:, -1], coeff)
            cocycles[dim].append(ccl)
    ret = {'dgms': dgms, 'cocycles': cocycles,
           'num_edges': res['num_edges'], 'dm': dm}
    return ret


def plot_dgms(diagrams, plot_only=None, title=None, xy_range=None,
              labels=None, colormap='default', size=20,
              ax_color=np.array([0.0, 0.0, 0.0]), colors=None,
              diagonal=True, lifetime=False, legend=True, show=False):
    """A helper function to plot persistence diagrams. 

    Parameters
    ----------

    diagrams: ndarray (n_pairs, 2) or list of diagrams
        A diagram or list of diagrams. If diagram is a list of diagrams, 
        then plot all on the same plot using different colors.

    plot_only: list of numeric
        If specified, an array of only the diagrams that should be plotted.

    title: string, default is None
        If title is defined, add it as title of the plot.

    xy_range: list of numeric [xmin, xmax, ymin, ymax]
        User provided range of axes. This is useful for comparing 
        multiple persistence diagrams.

    labels: string or list of strings
        Legend labels for each diagram. 
        If none are specified, we use H_0, H_1, H_2,... by default.

    colormap: string, default is 'default'
        Any of matplotlib color palettes. 
        Some options are 'default', 'seaborn', 'sequential'. 
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
        Plot life time of each point instead of birth and death. 
        Essentially, visualize (x, y-x).

    legend: bool, default is True
        If true, show the legend.

    show: bool, default is False
        Call plt.show() after plotting. If you are using self.plot() as part 
        of a subplot, set show=False and call plt.show() only once at the end.

    """

    if labels is None:
        # Provide default labels for diagrams if using self.dgm_
        labels = ["$H_0$", "$H_1$", "$H_2$", "$H_3$",
                  "$H_4$", "$H_5$", "$H_6$", "$H_7$", "$H_8$"]

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

    # Construct copy with proper type of each diagram
    # so we can freely edit them.
    diagrams = [dgm.astype(np.float32, copy=True) for dgm in diagrams]

    # find min and max of all visible diagrams
    concat_dgms = np.concatenate(diagrams).flatten()
    has_inf = np.any(np.isinf(concat_dgms))
    finite_dgms = concat_dgms[np.isfinite(concat_dgms)]

    if not xy_range:
        # define bounds of diagram
        ax_min, ax_max = np.min(finite_dgms), np.max(finite_dgms)
        ax_range = ax_max - ax_min

        # Give plot a nice buffer on all sides.
        # ax_range=0 when only one point,
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


class Rips(TransformerMixin):
    """sklearn class wrapper around Uli Bauer's ripser code 

    Parameters
    ----------
    maxdim : int, optional, default 1
        Maximum homology dimension computed. Will compute all dimensions 
        lower than and equal to this value. 
        For 1, H_0 and H_1 will be computed.

    thresh : float, default infinity
        Maximum distances considered when constructing filtration. 
        If infinity, compute the entire filtration.

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

    def __init__(self, maxdim=1, thresh=np.inf, coeff=2,
                 do_cocycles=False, verbose=True):
        self.maxdim = maxdim
        self.thresh = thresh
        self.coeff = coeff
        self.do_cocycles = do_cocycles
        self.verbose = verbose

        # Internal variables
        self.dgms_ = None
        self.cocycles_ = None
        self.dm_ = None  # Distance matrix
        self.metric_ = None
        self.num_edges_ = None  # Number of edges added

        if self.verbose:
            print("Rips(maxdim={}, thresh={}, coeff={}, do_cocycles={}, verbose={})".format(
                maxdim, thresh, coeff, do_cocycles, verbose))

    def transform(self, X, distance_matrix=False, metric='euclidean'):
        result = ripser(X, maxdim=self.maxdim, thresh=self.thresh,
                        coeff=self.coeff, do_cocycles=self.do_cocycles,
                        distance_matrix=distance_matrix, metric=metric)
        self.dgms_ = result['dgms']
        self.num_edges_ = result['num_edges']
        self.dm_ = result['dm']
        self.cocycles_ = result['cocycles']
        return self.dgms_

    def fit_transform(self, X, distance_matrix=False, metric='euclidean'):
        """
        Compute persistence diagrams for X data array and return the diagrams.

        Parameters
        ----------
        X: ndarray (n_samples, n_features)
            A numpy array of either data or distance matrix.

        distance_matrix: bool
            Indicator that X is a distance matrix, if not we compute a 
            distance matrix from X using the chosen metric.

        metric: string or callable
            The metric to use when calculating distance between instances in a 
            feature array. If metric is a string, it must be one of the options 
            specified in PAIRED_DISTANCES, including "euclidean", "manhattan", 
            or "cosine". Alternatively, if metric is a callable function, it is 
            called on each pair of instances (rows) and the resulting value 
            recorded. The callable should take two arrays from X as input and 
            return a value indicating the distance between them.

        Return
        ------
        dgms: list (size maxdim) of ndarray (n_pairs, 2)
            A list of persistence diagrams, one for each dimension less 
            than maxdim. Each diagram is an ndarray of size (n_pairs, 2) with 
            the first column representing the birth time and the second column 
            representing the death time of each pair.
        """
        self.transform(X, distance_matrix, metric)
        return self.dgms_

    def plot(self, diagrams=None, plot_only=None, title=None, xy_range=None,
             labels=None, colormap='default', size=20,
             ax_color=np.array([0.0, 0.0, 0.0]), colors=None, diagonal=True,
             lifetime=False, legend=True, show=True):
        """A helper function to plot persistence diagrams. 

        Parameters
        ----------
        diagrams: ndarray (n_pairs, 2) or list of diagrams
            A diagram or list of diagrams as returned from self.fit. 
            If diagram is None, we use self.dgm_ for plotting. 
            If diagram is a list of diagrams, then plot all on the same plot 
            using different colors.

        plot_only: list of numeric
            If specified, an array of only the diagrams that should be plotted.

        title: string, default is None
            If title is defined, add it as title of the plot.

        xy_range: list of numeric [xmin, xmax, ymin, ymax]
            User provided range of axes. This is useful for comparing 
            multiple persistence diagrams.

        labels: string or list of strings
            Legend labels for each diagram. 
            If none are specified, we use H_0, H_1, H_2,... by default.

        colormap: string, default is 'default'
            Any of matplotlib color palettes. 
            Some options are 'default', 'seaborn', 'sequential'. 
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
            Plot life time of each point instead of birth and death. 
            Essentially, visualize (x, y-x).

        legend: bool, default is True
            If true, show the legend.

        show: bool, default is True
            Call plt.show() after plotting. 
            If you are using self.plot() as part of a subplot, 
            set show=False and call plt.show() only once at the end.
        """

        if diagrams is None:
            # Allow using transformed diagrams as default
            diagrams = self.dgms_
        plot_dgms(diagrams, plot_only=plot_only, title=title,
                  xy_range=xy_range, labels=labels, colormap=colormap,
                  size=size, ax_color=ax_color, colors=colors,
                  diagonal=diagonal, lifetime=lifetime,
                  legend=legend, show=show)
