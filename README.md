# Ripser


Ripser is now a Python class implemented in the Scikit-learn style. It is easy to install, only requires that you have Cython installed first. It is even easier to use.

For the C++ library, see [Ripser/ripser](https://github.com/Ripser/ripser/releases/latest).

Details from the old readme can be found [here](docs/README.md).

## Setup

Installation requires Cython, and currently must be installed from source. An example of how to install is
```
git clone https://github.com/sauln/ripser
cd ripser
python -m venv venv
source venv/bin/activate
pip install Cython
pip install -e .
```

We use matplotlib for generating persistence diagrams


## Usage

```
import numpy as np
from ripser import Rips
r = Rips()

data = np.random.random((100,2))
diagram = r.fit_transform(data)
r.plot(diagram)
```

