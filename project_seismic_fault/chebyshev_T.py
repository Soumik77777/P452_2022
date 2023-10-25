import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.special

## Define Limits and Function
a, b = 1, 3

def input_func(y):
    return np.log(y)


def modi_func(y):
    y2 = ((a+b)/2) + ((b-a)*y/2)

    return np.sqrt(1-y**2) * input_func(y2) * (b-a)/2


deg = 100

cheb_return = np.polynomial.chebyshev.chebgauss(deg)

x = cheb_return[0]
w = cheb_return[1]

sum = 0
for i in range(len(x)):
    sum += w[i]*modi_func(x[i])

#print("Integrated by using Gauss-Chebyshev polynomials", sum)

print(scipy.special.expi())

