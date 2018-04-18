import pytest
import numpy as np

from ripser import Rips


class TestLibrary():
    # Does the library install in scope? Are the objects in scope?
    def test_import(self):
        import ripser
        from ripser import Rips
        assert 1

    def test_instantiate(self):
        rip = Rips()
        assert rip is not None

class TestParams():
    def test_defaults(self):
        data = np.random.random((100,3))

        rips = Rips()
        dgm = rips.fit_transform(data)

        assert len(dgm) == 2
        assert rips.coeff == 2


    def test_coeff(self):

        data = np.random.random((100,3))

        rips3 = Rips(coeff=3)
        dgm3 = rips3.fit_transform(data)

        rips2 = Rips(coeff=2)
        dgm2 = rips2.fit_transform(data)
        assert dgm2 is not dgm3, "This is a vacuous assertion, we only care that the above operations did not throw errors"
    
    def test_maxdim(self):
        data = np.random.random((100,3))

        # maxdim refers to the max H_p class, generate all less than

        rips0 = Rips(maxdim=0)
        dgm0 = rips0.fit_transform(data)
        assert len(dgm0) == 1

        rips1 = Rips(maxdim=1)
        dgm1 = rips1.fit_transform(data)
        assert len(dgm1) == 2

        rips2 = Rips(maxdim=2)
        dgm2 = rips2.fit_transform(data)
        assert len(dgm2) == 3


    def test_thresh(self):
        data = np.random.random((100,3))

        rips0 = Rips(thresh=0.1)
        dgm0 = rips0.fit_transform(data)

        rips1 = Rips(thresh=1)
        dgm1 = rips1.fit_transform(data)

        # Barcode of H_1 diagram will be smaller, right?
        assert len(dgm0[1]) < len(dgm1[1]), "Usually"