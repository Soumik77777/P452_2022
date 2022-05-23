import numpy as np
import library


def func(x):
    return np.exp(-x**2)

a = 653
m = 102809
n = 100000

seed = 77

without_samp = library.MonteCarlos(func, 0, 1, n, a, m, seed)
with_samp = library.montecarlo_imp(a, m, n, seed)


print("Integration value without impostance sampling= "+str(without_samp))
print("Integration value with importance sampling= "+str(with_samp))



'''
Integration value without impostance sampling= 0.7459327164223801
Integration value with importance sampling= 0.7455433612206332
'''


