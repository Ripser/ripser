
"""
To install the python bindings (on a mac)

    > make python
    > pip install -e .

"""


import pytest
import ripser


def test_import():
    import ripser

def test_instantiate():
    rip = ripser.Rips()

    assert rip is not None