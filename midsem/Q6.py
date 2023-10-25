###########

## Soumik Bhattacharyya, Roll- 1811155
## Q6. jacobi and gauss-seidel

import os
import csv
import math
import matplotlib.pyplot as plt


###########
def diff(xk0, xk1):
    sum = 0
    for i in range(len(xk0)):
        for j in range(len(xk0[i])):
            diff = abs(xk1[i][j] - xk0[i][j])
            sum = sum + diff
    return sum


def jacobi(A, B, eps, xk0):
    itr = 0
    xk1 = [[0] for x in range(len(A))]
    while diff(xk1, xk0) >= eps:
        if itr != 0:
            for i in range(len(xk1)):
                for j in range(len(xk1[i])):
                    xk0[i][j] = xk1[i][j]
        for i in range(len(A)):
            sum = 0
            for j in range(len(A[i])):
                if j != i:
                    sum = sum + (A[i][j] * xk0[j][0])
            xk1[i][0] = (1 / A[i][i]) * (B[i][0] - sum)
        itr = itr + 1
    return {'X': xk1, 'iterations': itr}


def gauss_seidel(A, B, eps, xk0):
    itr = 0
    xk1 = [[0] for x in range(len(A))]

    while diff(xk1, xk0) >= eps:
        if itr != 0:
            for i in range(len(xk1)):
                for j in range(len(xk1[i])):
                    xk0[i][j] = xk1[i][j]
        for i in range(len(A)):
            sum1 = 0
            sum2 = 0
            for j in range(i + 1, len(A[i])):
                sum2 = sum2 + (A[i][j] * xk0[j][0])
            for j in range(0, i):
                sum1 = sum1 + (A[i][j] * xk1[j][0])
            xk1[i][0] = (1 / A[i][i]) * (B[i][0] - sum1 - sum2)
        itr = itr + 1

    return {'X': xk1, 'iterations': itr}


matrix = [[-2,    0,    0,    -1,    0,    0.5],
		 [0,    4,    0.5,   0,    1,    0],
		 [0,    0.5,  1.5,   0,    0,    0],
		 [-1,    0,    0,    -2,    0,    1],
		 [0,    1,    0,     0,   -2.5,  0],
		 [0.5,  0,    0,     1,    0,   -3.75]]

vector = [[-1],
		 [0],
		 [2.75],
		 [2.5],
		 [-3],
		 [2]]

jacobi_soln = jacobi(matrix, vector, 0.00001, [[2], [3], [1], [1], [4], [2]])['X']
gs_soln = gauss_seidel(matrix, vector, 0.00001, [[2], [3], [1], [1], [4], [2]])['X']

print("The solution using Jacobi Method:")
print(jacobi_soln)

print("The solution using Gauss-seidel Method:")
print(gs_soln)






