.. ripser documentation master file, created by
   sphinx-quickstart on Sun Jul 22 20:37:23 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


|PyPI version| |Travis-CI| |Appveyor| |Codecov| |License: LGPL v3|

Ripser
========

This package provides the awesome Ripser project as an easy to use Python module. It is easy to install and is even easier to use.

For the original C++ library, see `Ripser/ripser <https://github.com/Ripser/ripser/releases/latest>`_.


Setup
------

Installation requires Cython, and currently must be installed from source. An example of how to install is

.. code:: python

    pip install Cython
    pip install Ripser


We use matplotlib for generating persistence diagrams


Usage
------

.. code:: python

    import numpy as np
    from ripser import ripser, plot_dgms

    data = np.random.random((100,2))
    diagrams = ripser(data)['dgms']
    plot_dgms(diagrams, show=True)



Note that there is also a <i>Rips</i> object with the same functionality, which conforms to the Scikit-learn API.

.. code:: python

    import numpy as np
    from ripser import Rips
    r = Rips()

    data = np.random.random((100,2))
    diagram = r.fit_transform(data)
    r.plot(diagram, show=True)



.. |PyPI version| image:: https://badge.fury.io/py/ripser.svg
   :target: https://badge.fury.io/py/ripser

.. |Travis-CI| image:: https://travis-ci.org/ctralie/ripser.svg?branch=master
    :target: https://travis-ci.org/ctralie/ripser

.. |Appveyor| image:: https://ci.appveyor.com/api/projects/status/sfy7yybs66e5qanu?svg=true
    :target: https://ci.appveyor.com/project/ctralie/ripser
.. |Codecov| image:: https://codecov.io/gh/ctralie/ripser/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/ctralie/ripser
.. |License: LGPL v3| image:: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/lgpl-3.0



.. toctree::
    :maxdepth: 2
    :caption: Contents:

    Ripser Demonstration
    Approximate Sparse Filtrations
    Sparse Distance Matrices
    Representative Cocycles
    Lower Star Time Series


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
