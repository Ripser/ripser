""" Moebius strip example.  A fake H1 class exists with Z2 coefficients
    which does not exist with Z3 coefficients
"""

import numpy as np
from ripser import Rips

t = np.linspace(0, 2*np.pi, 500)
X = np.zeros((500,3))
X[:,0] = (2 + np.cos(t)) * np.cos(2*t)
X[:,1] = (2 + np.cos(t)) * np.sin(2*t)
X[:,2] = np.sin(t)


print('Computing over Z2...')
r = Rips(maxdim=2, coeff=2)
dgmZ2 = r.fit_transform(X)
print('Finished Z2')

print('Computing over Z3...')
r = Rips(maxdim=2, coeff=2)
dgmZ3 = r.fit_transform(X)
print('Finished Z3')

r.plot(dgmZ2)
r.plot(dgmZ3)
