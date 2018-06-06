[![PyPI version](https://badge.fury.io/py/ripser.svg)](https://badge.fury.io/py/ripser)
[![Build Status](https://travis-ci.org/ctralie/ripser.svg?branch=master)](https://travis-ci.org/ctralie/ripser)
[![codecov](https://codecov.io/gh/ctralie/ripser/branch/master/graph/badge.svg)](https://codecov.io/gh/ctralie/ripser)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

# Ripser


Ripser is now a Python class. It is easy to install, only requires that you have Cython installed first. It is even easier to use.

For the C++ library, see [Ripser/ripser](https://github.com/Ripser/ripser/releases/latest).

Details from the old readme can be found [here](docs/README.md).

## Setup

Installation requires Cython, and currently must be installed from source. An example of how to install is
```
pip install Cython
pip install Ripser
```

We use matplotlib for generating persistence diagrams


## Usage

```
import numpy as np
from ripser import ripser, plot_dgms

data = np.random.random((100,2))
diagrams = ripser(data)['dgms']
plot_dgms(diagrams, show=True)
```


Note that there is also a <i>Rips</i> object with the same functionality, which conforms to the Scikit-learn style

```
import numpy as np
from ripser import Rips
r = Rips()

data = np.random.random((100,2))
diagram = r.fit_transform(data)
r.plot(diagram, show=True)
```

