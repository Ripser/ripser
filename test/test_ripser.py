import pytest
import numpy as np

from ripser import ripser
from sklearn import datasets
from sklearn.metrics.pairwise import pairwise_distances
from scipy import sparse


def makeSparseDM(X, thresh):
    """
    Helper function to make a sparse distance matrix
    """
    N = X.shape[0]
    D = pairwise_distances(X, metric='euclidean')
    [I, J] = np.meshgrid(np.arange(N), np.arange(N))
    I = I[D <= thresh]
    J = J[D <= thresh]
    V = D[D <= thresh]
    return sparse.coo_matrix((V, (I, J)), shape=(N, N)).tocsr()


class TestLibrary():
    # Does the library install in scope? Are the objects in scope?
    def test_import(self):
        import ripser
        from ripser import ripser, plot_dgms
        assert 1


class TestTransform():
    def test_input_warnings(self):
        data = np.random.random((3, 10))

        with pytest.warns(UserWarning, match='has more columns than rows') as w:
            ripser(data)

        data = np.random.random((3, 3))
        with pytest.warns(UserWarning, match='input matrix is square, but the distance_matrix') as w:
            ripser(data)

    def test_non_square_dist_matrix(self):
        data = np.random.random((3, 10))

        with pytest.raises(Exception):
            ripser(data, distance_matrix=True)


class TestParams():
    def test_defaults(self):
        data = np.random.random((100, 3))
        dgms = ripser(data)['dgms']
        assert len(dgms) == 2

    def test_coeff(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        dgm3 = ripser(data, coeff=3)['dgms']
        dgm2 = ripser(data)['dgms']
        assert dgm2 is not dgm3, "This is a vacuous assertion, we only care that the above operations did not throw errors"

    def test_maxdim(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        # maxdim refers to the max H_p class, generate all less than
        dgms0 = ripser(data, maxdim=0)['dgms']
        assert len(dgms0) == 1

        dgms1 = ripser(data)['dgms']
        assert len(dgms1) == 2

        dgms2 = ripser(data, maxdim=2)['dgms']
        assert len(dgms2) == 3

    def test_thresh(self):
        np.random.seed(3100)
        data = np.random.random((100, 3))

        dgms0 = ripser(data, thresh=0.1)['dgms']
        dgms1 = ripser(data)['dgms']

        # Barcode of H_1 diagram will be smaller, right?
        assert len(dgms0[1]) < len(dgms1[1]), "Usually"

    def test_sparse(self):
        np.random.seed(10)
        thresh = 1.1

        # Do dense filtration with threshold
        data = datasets.make_circles(n_samples=100)[
            0] + 5 * datasets.make_circles(n_samples=100)[0]
        res0 = ripser(data, thresh=thresh)

        # Convert to sparse matrix first based on threshold,
        # then do full filtration
        D = makeSparseDM(data, thresh)
        res1 = ripser(D, distance_matrix=True)

        # The same number of edges should have been added
        assert res0['num_edges'] == res1['num_edges']

        dgms0 = res0['dgms']
        dgms1 = res1['dgms']
        I10 = dgms0[1]
        I11 = dgms1[1]
        idx = np.argsort(I10[:, 0])
        I10 = I10[idx, :]
        idx = np.argsort(I11[:, 0])
        I11 = I11[idx, :]
        assert np.allclose(I10, I11)
