import math
import csv
import random
import numpy as np
import scipy
import scipy.special


######   For integrating/ decomposing using quadrature

## Type 1 / T

def modi_func_T(input_func, y, a, b):             ## func is multiplied by w, limit is changed
    y2 = ((a+b)/2) + ((b-a)*y/2)
    return input_func(y2) * ((b-a)/2) * np.sqrt(1-y**2)

def chebyshev_T_sw(n):                          ## returns s and w lists of a func in a given limit and degree
    s, w = [], []
    for i in range(1,n+1):
        s.append(np.cos( (np.pi * (i-0.5)) / n))
        w.append(np.pi / n)
    return s, w


def chebyshev_T_complimentary(n):
    r = []
    for i in range(1, n):
        r.append(np.cos(np.pi * i / n))
    return r




###########

## Type 2 / U

def modi_func_U(input_func, y, a, b):             ## func is multiplied by w, limit is changed
    y2 = ((a+b)/2) + ((b-a)*y/2)
    return input_func(y2) * ((b-a)/2) / np.sqrt(1-y**2)


def chebyshev_U_sw(n):                          ## returns s and w lists of a func in a given limit and degree
    s, w = [], []
    for i in range(1,n+1):
        k = np.pi / (n + 1)
        s.append(np.cos(i * k))
        w.append(k * (np.sin(i * k) ** 2))
    return s, w


def chebyshev_U_complimentary(n):
    r = []
    for i in range(1, n+2):
        k = np.pi * (i - 0.5) / (n+1)
        r.append(np.cos(k))
    return r



def cheb_type2(func, type=T):                           ## returns integrated value using func, limits, degree
    print("Please check the input function.")
    def check_func():
        print("Type 'yes' to proceed. Type 'no' to break the loop.")
        func_check = input(str(": "))
        return func_check

    if check_func() == 'yes':
        l = float(input("Lower bound of integration: "))
        u = float(input("Upper bound of integration: "))
        deg = int(input("Degree of Chebyshev Polynomial: "))

        s, w = chebyshev_U_sw(deg)

        sum = 0
        for i in range(len(s)):
            sum += w[i] * modi_func_U(func, s[i], l, u)

        return sum


    elif check_func() == 'no':
        print("Modify the function and rerun")
        return 0
    else:
        check_func()




########  Functions in problem



def func_E(r, s):
    U = (r + 1) / (s + 1)
    k = 2 * np.sqrt(U) / (1 + U)
    E = scipy.special.ellipe(k)
    return E

def func_F(r,s):
    U = (r + 1) / (s + 1)
    k = 2 * np.sqrt(U) / (1 + U)
    F = scipy.special.ellipk(k)
    return F

def coeff_slip(r, s):
    term1 = func_E(r,s) / (s-r)
    term2 = func_F(r,s) / (s + r + 2)

    return term1 + term2


def lambda_func(lam, r):
    term = (lam**2) * ((r + 1) **2)
    return scipy.special.expi(term)


def jaco_lambda(lam, r):
    return 2 * np.exp(-((lam**2) * ((r+1) ** 2))) / lam


















################   Newton-Raphson

def derivative (f, x, h):
    g = (f(x+h) - f(x-h))/(2*h)
    return g

def NewtonRaphson(f, x, h, e):
    g = derivative(f, x, h)
    d = f(x) / g
    list_x = [x]

    while abs(d) >= e:
        g = derivative(f, x, h)
        d = f(x) / g
        x = x - d
        list_x.append(x)

    return (x, list_x)

############## Matrix




