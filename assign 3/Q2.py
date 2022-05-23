import matplotlib.pyplot as plt
import numpy as np

import library


def f(x, y, z):
    return -n ** 2 * np.pi ** 2 * y

n = 1
X1, Y1 = library.shooting(f, 0, 1, 0, 0)
n = 2
X2, Y2 = library.shooting(f, 0, 1, 0, 0)

Y1 = library.normalize(Y1)
Y2 = library.normalize(Y2)

plt.plot(X1, Y1, label='Lowest Energy State')
plt.plot(X2, Y2, color='red', label='Second Lowest Energy State')
plt.axvline(x=0, color='green', lw=0.7)
plt.axvline(x=1, color='green', lw=0.7)
plt.xlabel('x')
plt.ylabel('$\psi (x)$')
plt.legend()
plt.grid()
plt.show()