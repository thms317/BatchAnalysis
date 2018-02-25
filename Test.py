# from Definitions import WLC
from definitions import *
import numpy as np
import matplotlib.pyplot as plt

L = 3646  # contour length (bp)
p = 50  # persistence length (nm)
S = 1000  # stretch modulus (pN)
x0 = 0  # offset (nm)

force = []
extension = []

for f in range(1, 500):
    f /= 10
    force.append(f)
    extension.append(WLC(f, p, L, S, x0))

plt.plot(extension, force)
plt.show()
