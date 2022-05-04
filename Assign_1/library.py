import csv
import os
import copy
import math
import random
from fractions import Fraction
import numpy as np




#### Gauss Jordan -------------------------------
def partial_pivot (m, v):
    n = len(v)

    for i in range (n-1):
        if m[i][i] ==0:
            for j in range (i+1,n):
                if abs(m[j][i]) > abs(m[i][i]):
                    m[i], m[j] = m[j], m[i]
                    v[i], v[j] = v[j], v[i]
    return (m,v)


def gauss_jordan(mat_M, b):
    n = len(mat_M)
    #do partial pivoting
    partial_pivot(mat_M, b)

    for r in range(n):
        #make the diagonal element 1
        pivot = mat_M[r][r]
        for c in range(r,n):
            mat_M[r][c] = mat_M[r][c]/pivot
        b[r] = b[r]/pivot

        #make the other element in that column 0
        for r1 in range(n):
            #nothing to do for the diagonal element or if it already is 0
            if (r1 == r) or (mat_M[r1][r] == 0):
                continue
            else:
                factor = mat_M[r1][r]
                for c in range(r,n):
                    mat_M[r1][c] = mat_M[r1][c] - factor*mat_M[r][c]
                b[r1] = b[r1] - factor*b[r]
    return b



## LU Decompose ---------------------------
def ludecompose (matrix, n):

    upper_mat = [[0 for i in range(n)] for j in range(n)]
    lower_mat = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):

        for j in range(i, n): #calculating upper matrix
            sum = 0
            for k in range(i):
                sum += (lower_mat[i][k] * upper_mat[k][j])
            upper_mat[i][j] = matrix[i][j] - sum

        for j in range(i, n): #calculating lower matrix
            if (i == j):
                lower_mat[i][i] = 1
            else:
                sum = 0
                for k in range(i):
                    sum += (lower_mat[j][k] * upper_mat[k][i])

                lower_mat[j][i] = ((matrix[j][i] - sum) / upper_mat[i][i])

    return (lower_mat, upper_mat)


def forward_backward_substitution (lower_mat, upper_mat, vector, n):
    '''
    If we have LUx=B,
    first we solve Ly=B, then Ux=y
    '''
    # forward-substitution
    y = [0] * n
    for i in range(n):
        sum = 0
        for j in range(i):
            sum += lower_mat[i][j] * y[j]

        y[i] = vector[i] - sum

    #backward-substitution
    x = [0] * n
    for i in reversed(range(n)):
        sum = 0
        for j in range(i + 1, n):
            sum+= upper_mat[i][j] * x[j]
        x[i] = (y[i] - sum)/ upper_mat[i][i]

    return (x)

def LUSolve(mat,vec):
    n = len(mat)
    L, U = ludecompose(mat, n)
    x = forward_backward_substitution(L,U,vec,n)
    return x



## Gauss Jacobi -------------------------------

def gauss_jacobi(a, b):
    m = len(b[0])
    n = len(a)
    epsilon = math.pow(10, -4)

    x_k = [[1 for y in range(m)] for x in range(n)]
    x_k1 = [[0 for y in range(m)] for x in range(n)]
    norm = 1
    while norm > epsilon:
        norm = 0
        for i in range(n):
            for j in range(m):
                inner_sum = 0
                for k in range(n):
                    if k != i:
                        inner_sum = inner_sum + a[i][k] * x_k[k][j]
                x_k1[i][j] = (1 / a[i][i]) * (b[i][j] - inner_sum)
            for j in range(m):
                norm += math.pow((x_k1[i][j] - x_k[i][j]), 2)
        norm = math.pow(norm, 0.5)
        x_k = copy.deepcopy(x_k1)
    return x_k1


def multiply(ain,b) :
	m=len(b)
	k2 = [0 for y in range(m)]
	for i in range (m) :
		for j in range (m) :
			k2[i] = k2[i] + (ain[i][j]*b[j])
	return k2





## Gauss Seidel -----------------------------------

def gauss_siedel(a, b):
    m = len(b[0])
    n = len(a)
    epsilon = math.pow(10, -5)
    x_k = [[1 for y in range(m)] for x in range(n)]
    norm = 1
    while norm > epsilon:
        norm = 0
        for i in range(n):
            for j in range(m):
                inner_sum = 0
                for k in range(n):
                    if k != i:
                        inner_sum = inner_sum + a[i][k] * x_k[k][j]
                l = x_k[i][j]
                x_k[i][j] = (1 / a[i][i]) * (b[i][j] - inner_sum)
                norm += math.pow((x_k[i][j] - l), 2)
        norm = math.pow(norm, 0.5)
    return x_k




## Conjugate Gradient ----------------------------

def conjGrad(A,b,tol):
    xk = []
    for i in range(len(A)):
        xk.append(0)
    rk = np.dot(A, xk) - b
    pk = -rk
    rk_norm = np.linalg.norm(rk)    
    n = 0    
    while rk_norm > tol:
        Apk = np.dot(A, pk)
        rkrk = np.dot(rk, rk)        
        alpha = rkrk / np.dot(pk, Apk)
        xk = xk + alpha * pk
        rk = rk + alpha * Apk
        beta = np.dot(rk, rk) / rkrk
        pk = -rk + beta * pk        
        n = n+1        
        rk_norm = np.linalg.norm(rk)         
    return xk

def invCG(A,tol,identity):
    n = len(A)
    I = identity
    Inv = [[0 for i in range(n)]for j in range(n)]
    cb =  [[0 for i in range(1)]for j in range(n)]
    ci =  [[0 for i in range(1)]for j in range(n)]
    for i in range(n):
        for j in range(len(I)):
            ci[j][0] = I[j][i]
        c = np.array(ci)
        r = c.T
        f = r[0].tolist()
        b = conjGrad(A,f,tol)
        for k in range(n):
            Inv[k][i] = b[k]
    return Inv
