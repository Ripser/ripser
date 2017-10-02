import ripser
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 2*np.pi, 100)
x = np.cos(t)
D = np.abs(x[:, None] - x[None, :])

modulus = 3
dim_max = 1
threshold = np.max(D)*2

print(ripser.ripser(D, modulus, dim_max, threshold))
