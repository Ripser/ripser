""" Moebius strip example.  A fake H1 class exists with Z2 coefficients
    which does not exist with Z3 coefficients
"""

import numpy as np
from ripser import Rips
import matplotlib.pyplot as plt

N = 50
t = np.linspace(0, 2*np.pi, N)
X = np.zeros((N,3))
X[:,0] = (2 + np.cos(t)) * np.cos(2*t)
X[:,1] = (2 + np.cos(t)) * np.sin(2*t)
X[:,2] = np.sin(t)


print('Computing over Z2...')
r = Rips(maxdim=1, coeff=2)
dgmZ2 = r.fit_transform(X)
print('Finished Z2')

print('Computing over Z3...')
r = Rips(maxdim=1, coeff=3)
dgmZ3 = r.fit_transform(X)
print('Finished Z3')


# Plot the two diagrams
plt.figure(1)

# linear
plt.subplot(121)
r.plot(dgmZ2, show=False)
plt.title("Diagram over Z2")

plt.subplot(122)
r.plot(dgmZ3, show=False)
plt.title("Diagram over Z3")
plt.show()
