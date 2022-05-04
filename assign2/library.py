import copy
import math
import random
from fractions import Fraction
import matplotlib.pyplot as plt
import numpy as np


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

"""---------------------------------------------------------------------------------------"""
def constructM(con, r, c):
    A=[]
    B=[]
    with open(con,'r+') as file:
        i=0
        while i<r:
            a=file.readline()
            b=a.split()
            j=0
            while j<c:
                B.append(float(b[j]))
                j=j+1
            A.append(B)
            B=[]
            i=i+1
        return A



# Reading a Matrix from a file
def read_Mat(fil, A):
    file = open(fil, 'r')
    for line in file:
        ns = line.split()
        no = [float(Fraction(n)) for n in ns]
        A.append(no)
    file.close()

def makmat(file,r,c):
    A=[]
    B=[]
    with open(file,'r+') as file:
        i=0
        while i<r:
            x=file.readline().strip()

            y=x.split(' ')
            j=0
            while j<c:
                B.append(float(y[j]))
                j=j+1
            A.append(B)
            B=[]
            i=i+1
        return A


def read_MatCol(fil):
    file = open(fil, 'r')
    x = []
    y = []
    std = []
    for line in file:
        ns = line.split()
        no = [float(n) for n in ns]
        x += [[no[0]]]
        y += [[no[1]]]
        std += [[no[2]]]
    file.close()
    return x, y, std


"""---------------------------------------------------------------------------------------"""


# To print the Matrix
def print_Mat(x):
    for r in range(len(x)):
        print(x[r])


"""---------------------------------------------------------------------------------------"""


# Factorial Method
def factorial(n):
    fact = 1
    if n < 0:
        print("Sorry, factorial does not exist for negative numbers")
    elif n == 0:
        return 1
    else:
        for i in range(1, n + 1):
            fact = fact * i
        return fact


"""---------------------------------------------------------------------------------------"""


# Function for partial pivoting the Augmented Matrix / Only Matrix
def partial_pivot(a, b):
    n = len(a)
    counter = 0
    for r in range(0, n):
        # if abs(a[r][r]) == 0:
        for r1 in range(r + 1, n):
            if abs(a[r1][r]) > abs(a[r][r]):
                counter = counter + 1
                for x in range(0, n):
                    d1 = a[r][x]
                    a[r][x] = a[r1][x]
                    a[r1][x] = d1
                if b != 0:
                    d1 = b[r]
                    b[r] = b[r1]
                    b[r1] = d1
    return counter


"""---------------------------------------------------------------------------------------"""


def mat_Add(a, b):
    m = len(b[0])
    n = len(a)
    rmatrix = [[0 for y in range(m)] for x in range(n)]
    for i in range(len(a)):
        for j in range(len(b[0])):
            rmatrix[i][j] = a[i][j] + b[i][j]
    return rmatrix


def mat_Subtract(a, b):
    m = len(b[0])
    n = len(a)
    rmatrix = [[0 for y in range(m)] for x in range(n)]
    for i in range(len(a)):
        for j in range(len(b[0])):
            rmatrix[i][j] = a[i][j] - b[i][j]
    return rmatrix


# Multiplication of Matrices
def mat_Multiply(a, b):
    m = len(b[0])
    l = len(b)
    n = len(a)
    p2 = [[0 for y in range(m)] for x in range(n)]
    for x in range(n):
        for i in range(m):
            for y in range(l):
                p2[x][i] = p2[x][i] + (a[x][y] * b[y][i])
    return p2


#new addition

"""---------------------------------------------------------------------------------------"""


def mat_transpose(m):
    transpose = [[m[j][i] for j in range(len(m))] for i in range(len(m[0]))]
    return transpose


"""---------------------------------------------------------------------------------------"""


# Gauss-Jordan Method
def gauss_Jordan(a, b):
    n = len(a)
    bn = len(b[0])
    for r in range(0, n):
        partial_pivot(a, b)
        pivot = a[r][r]
        for c in range(r, n):
            a[r][c] = a[r][c] / pivot
        for c in range(0, bn):
            b[r][c] = b[r][c] / pivot
        for r1 in range(0, n):
            if r1 == r or a[r1][r] == 0:
                continue
            else:
                factor = a[r1][r]
                for c in range(r, n):
                    a[r1][c] = a[r1][c] - factor * a[r][c]
                for c in range(0, bn):
                    b[r1][c] = b[r1][c] - factor * b[r][c]


"""---------------------------------------------------------------------------------------"""


# Forward- Backward Substitution
def forwardbackward_Substitution(a, b):
    m = len(b[0])
    n = len(a)
    # forward substitution
    y = [[0 for y in range(m)] for x in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0
            for k in range(i):
                s = s + a[i][k] * y[k][j]
            y[i][j] = b[i][j] - s
    # backward substitution
    x = [[0 for y in range(m)] for x in range(n)]
    for i in range(n - 1, -1, -1):
        for j in range(m):
            s = 0
            for k in range(i + 1, n):
                s = s + a[i][k] * x[k][j]
            x[i][j] = (y[i][j] - s) / a[i][i]

    return x


"""---------------------------------------------------------------------------------------"""


# L-U decomposition
def lu_Decomposition(a):
    n = len(a)
    for i in range(n):
        for j in range(n):
            s = 0
            if i <= j:
                for k in range(i):
                    s = s + (a[i][k] * a[k][j])
                a[i][j] = a[i][j] - s
            else:
                for k in range(j):
                    s = s + (a[i][k] * a[k][j])
                a[i][j] = (a[i][j] - s) / a[j][j]


"""---------------------------------------------------------------------------------------"""


# Finding determinant of Upper Triangular Matrix
def uppertriangular_Determinant(a):
    n = len(a)
    p = 1
    for i in range(n - 2):
        p = p * a[i][i]
    p = p * (a[n - 2][n - 2] * a[n - 1][n - 1])
    return p


# Finding first derivative
def first_Derivative(f, x):
    h = math.pow(10, -4)
    return (f(x + h) - f(x - h)) / (2 * h)


def polynomialsecond_Derivative(p, c, i, x):
    h = math.pow(10, -6)
    m = polynomialfirst_Derivative(p, c, i, x + h)
    n = polynomialfirst_Derivative(p, c, i, x - h)
    return (m - n) / (2 * h)

#mean of a array of numbers
def mean(a):
    n = len(a)
    s = 0
    for i in range(n):
        s = s + a[i]
    return s / n

#linear fit
def StraightLine_Fitting(x, y):
    n = len(x)
    xbar = mean(x)
    ybar = mean(y)
    Sxx, Syy, Sxy = 0, 0, 0
    for i in range(n):
        Sxx = Sxx + math.pow((x[i] - xbar), 2)
        Syy = Syy + math.pow((y[i] - ybar), 2)
        Sxy = Sxy + ((x[i] - xbar) * (y[i] - ybar))
    sigmax2 = Sxx / n
    sigmay2 = Syy / n
    covariance = Sxy / n
    b = covariance / sigmax2

    a = ybar - b * xbar

    #Correlation Coefficient
    r2 = (Sxy ** 2) / (Sxx * Syy)
    r = math.sqrt(r2)
    return a, b, r


##Simultaneous Linear Equations solver
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
"""---------------------------------------------------------------------------------------"""


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


"""---------------------------------------------------------------------------------------"""




"""---------------------------------------------------------------------------------------"""


# power method
def power_method(a, eigen_position):
    epsilon = math.pow(10, -4)
    n = len(a)
    x0 = [[1] for x in range(n)]
    x1 = matrix_Multiplication(a, x0)
    l1 = matrix_Multiplication(matrix_transpose(x1), x0)[0][0] / matrix_Multiplication(matrix_transpose(x0), x0)[0][0]
    x0 = copy.deepcopy(x1)
    l2 = l1
    condition = True
    while condition:
        x1 = matrix_Multiplication(a, x0)
        l2 = matrix_Multiplication(matrix_transpose(x1), x0)[0][0] / matrix_Multiplication(matrix_transpose(x0), x0)[0][
            0]
        if abs(l2 - l1) < epsilon:
            condition = False
        else:
            x0 = copy.deepcopy(x1)
            l1 = l2
    x1 = copy.deepcopy(x0)
    norm_x1 = math.sqrt(matrix_Multiplication(matrix_transpose(x1), x1)[0][0])
    x1 = [[j / norm_x1 for j in i] for i in x1]
    if eigen_position == 1:
        print("eigen value = ", l2)
        print("Normalized eigen vector = ", x1)
    elif eigen_position == 2:
        u = [[l2 * j for j in i] for i in matrix_Multiplication(x1, matrix_transpose(x1))]
        a = matrix_Subtraction(a, u)
        power_method(a, 1)


"""---------------------------------------------------------------------------------------"""
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

# functions for jacobi and givens rotation for eigenvalue
def max_off_diagonal(a):
    max_offdiagonal = a[0][1]
    k = 0
    l = 1
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if abs(a[i][j]) > abs(max_offdiagonal):
                max_offdiagonal = a[i][j]
                k = i
                l = j
    return max_offdiagonal, k, l


def jacobi_givens(a):
    epsilon = math.pow(10, -4)
    n = len(a)
    maxi, k, l = max_off_diagonal(a)
    while abs(maxi) >= epsilon:
        rotation = [[(1 if i == j else 0) for j in range(n)] for i in range(n)]
        if a[k][k] - a[l][l] == 0:
            theta = math.pi / 4
        else:
            theta = math.atan((2 * a[k][l]) / (a[l][l] - a[k][k])) / 2
        rotation[k][k] = math.cos(theta)
        rotation[l][l] = math.cos(theta)
        rotation[k][l] = math.sin(theta)
        rotation[l][k] = -1 * math.sin(theta)
        # a_t = matrix_Multiplication(matrix_transpose(rotation),
        # matrix_Multiplication(a_t, rotation))
        # print("theta = ", theta*(180/math.pi))
        temp0 = a[k][k]
        temp1 = a[l][l]
        a[k][k] = temp0 * (math.cos(theta) ** 2) + temp1 * (math.sin(theta) ** 2) - 2 * a[k][l] * math.cos(
            theta) * math.sin(theta)
        a[l][l] = temp1 * (math.cos(theta) ** 2) + temp0 * (math.sin(theta) ** 2) + 2 * a[l][k] * math.cos(
            theta) * math.sin(theta)
        a[k][l] = (temp0 - temp1) * math.cos(theta) * math.sin(theta) + a[k][l] * (
                math.cos(theta) ** 2 - math.sin(theta) ** 2)
        a[l][k] = a[k][l]
        for i in range(n):
            temp0 = a[i][k]
            temp1 = a[i][l]
            if i != k and i != l:
                a[i][k] = temp0 * math.cos(theta) - temp1 * math.sin(theta)
                a[k][i] = a[i][k]
                a[i][l] = temp1 * math.cos(theta) + temp0 * math.sin(theta)
                a[l][i] = a[i][l]
        maxi, k, l = max_off_diagonal(a)
    # print(a)
    print("Eigen Values are: ")
    [print(a[i][i]) for i in range(n)]


"""---------------------------------------------------------------------------------------"""


# Chi Square Fit
# ~ def to_log(a):
    # ~ y11 = [[0] for i in range(len(a))]
    # ~ for i in range(len(a)):
        # ~ y11[i][0] = math.log(a[i][0])
    # ~ return y11


# ~ def from_log(a):
    # ~ y11 = [[0] for i in range(len(a))]
    # ~ for i in range(len(a)):
        # ~ y11[i][0] = math.exp(a[i][0])


# ~ def chisquarefitting(x, y, std, degree_n):
    # ~ N = len(x)
    # ~ A = [[0 for j in range(degree_n + 1)] for i in range(degree_n + 1)]
    # ~ b = [[0] for i in range(degree_n + 1)]
    # ~ for i in range(N):
        # ~ for j in range(degree_n + 1):
            # ~ b[j][0] += (x[i][0] ** j) * y[i][0] / (std[i][0] ** 2)
            # ~ for k in range(degree_n + 1):
                # ~ A[j][k] += (x[i][0] ** (j + k)) / (std[i][0] ** 2)
    # ~ lu_Decomposition(A)
    # ~ coeff = forwardbackward_Substitution(A, b)
    # ~ # print(coeff)
    # ~ chi = 0
    # ~ for i in range(N):
        # ~ sum = 0
        # ~ for j in range(degree_n + 1):
            # ~ sum += coeff[j][0] * math.pow(x[i][0], j)
        # ~ chi += math.pow((y[i][0] - sum) / std[i][0], 2)
    # ~ print("chi_square/nu =", chi / (N - degree_n))
    # ~ I = [[(1 if i == j else 0) for j in range(degree_n + 1)] for i in range(degree_n + 1)]
    # ~ A_inverse = forwardbackward_Substitution(A, I)
    # ~ return coeff, A_inverse, chi


# ~ # Exponential Fit converts the required values to logarithmic scale and vice versa
# ~ def Exponential_chi_fit(x, y, std, degree_n):
    # ~ y1 = to_log(y)
    # ~ x = to_log(x)
    # ~ std11 = [[0] for i in range(len(std))]
    # ~ for i in range(len(std)):
        # ~ std11[i][0] = std[i][0] / y[i][0]
    # ~ coeff, A_inverse, chi = chi_square_fitting(x, y1, std11, degree_n)
    # ~ k = math.exp(coeff[0][0])
    # ~ print(coeff)
    # ~ print("coeff[0] = ", k)
    # ~ print("Errors in Coefficients =", k * math.sqrt(A_inverse[0][0]))
    # ~ for i in range(1, len(A_inverse)):
        # ~ print("Errors in Coefficients =", math.sqrt(A_inverse[i][i]))
def ChiSquareFit(x,y,sig):     
    n = len(sig)
    s, sx, sy, sxx, sxy, syy = 0,0,0,0,0,0
    for i in range(n):
        var = (1/sig[i]**2)
        s = s + (1/sig[i]**2)
        sx = sx + x[i]*var
        sy = sy + y[i]*var
        sxx = sxx + (x[i]**2)*var
        syy = syy + (y[i]**2)*var
        sxy = sxy + x[i]*y[i]*var
    Del = s*sxx - (sx)**2
    a = ((sxx*sy) - (sx*sxy))/Del       
    b = ((s*sxy) - (sx*sy))/Del         
    sigmaa = sxx/Del                        
    sigmab = s/Del
    covarianceb = -sx/Del
    r2 = sxy/(sxx*syy)
    dof = n - 2                             
    return a,b,sigmaa,sigmab,covarianceb,r2,dof

def Plot_Analytical(f, a, b, h):
    X = []
    Y = []
    x = a
    while x<b:
        X.append(x)
        Y.append(f(x))
        x = x+h
    return X,Y
    

# JackKnife Method

def jackknife_estimation(function, data):
    n = len(data)
    f_JK = [[0] for i in range(n)]
    f_mean_JK = 0
    for i in range(n):
        sum = 0
        for j in range(n):
            if i != j:
                sum += data[j][0]
        f_JK[i][0] = function(sum / (n - 1))
        f_mean_JK += f_JK[i][0]

    f_mean_JK = f_mean_JK / n

    # For standard error

    f_sigma_square = 0
    for i in range(n):
        f_sigma_square += math.pow(f_JK[i][0] - f_mean_JK, 2)
    f_sigma_square = f_sigma_square * (n - 1) / n

    print(f_mean_JK, f_sigma_square ** 0.5)


def powerMethodIndividual(A,tol):
    x = [0]*(len(A))
    x[0] = 1
    oldEigenVal = 1
    eigenVal = 0
    while abs(oldEigenVal-eigenVal)>tol:
            oldEigenVal=eigenVal
            x = np.dot(A,x)
            eigenVal = max(abs(x))
            x = x/eigenVal
    return eigenVal,x

def powerMethod(A,tol):
    n = len(A)
    eig = []
    E,V = powerMethodIndividual(A,tol)
    eig.append(E)
    if n>1:
        count = n-1
        while count != 0:
            V = V/np.linalg.norm(V)
            V = np.array([V])
            A = A - (E*V*V.T)
            E,V = powerMethodIndividual(A,tol)
            eig.append(E)
            count = count - 1 
    return eig



def pseudoRandoms(a,m,n,seed):    
    x = [seed]
    for i in range(n):
        temp = (a*x[i])%m
        x.append(temp)
    x = [i/m for i in x]
    return x


def MonteCarlos(f,a,b,N,aa,m,seed):
    h = (b-a)/float(N)
    I = 0     
    xk = pseudoRandoms(aa,m,N,seed)
    for i in range(N):
        X = a + (b-a)*xk[i]
        I = I + f(X)        #Caculate integral
    return h*I



def sum1(X, n):
    n = n + 1
    suMatrix = []
    j = 0
    while j<2*n:
        sums = 0
        i = 0
        while i< len(X):
            sums = sums + (X[i])**j
            i = i + 1
        suMatrix.append(sums)
        j = j+1
    return suMatrix

def makemat(suMatrix, n):
    n = n + 1
    m = [[0 for i in range(n)]for j in range(n)]
    i = 0
    while i<n:
        j = 0
        while j<n:
            m[i][j] = suMatrix[j+i]
            j = j+1
        i = i + 1
    return m
def sum2(X, Y, n):
    n = n+1
    suMatrix = []
    j = 0
    while j<n:
        sums = 0
        i = 0
        while i< len(X):
            sums = sums + ((X[i])**j)*Y[i]
            i = i + 1
        suMatrix.append(sums)
        j = j+1
    return suMatrix

def fit(X,Y,degree):
    k = sum1(X, degree)        
    m = makemat(k, degree)      
    qw = sum2(X, Y, degree)
    X = gauss_jordan(m, qw)
    return X[0], X[1], X[2], X[3]


def basis(deg,x):
    if (deg == 0):
        return (1)
    elif (deg == 1):
        return (2*x - (1))
    elif (deg == 2):
        return ((8*x**2) - (8*x) + 1)
    elif (deg == 3):
        return ((32*x**3) - (48*x**2) + (18*x) - 1)
    elif (deg == 4):
        return (2*x - (1))
    elif (deg == 5):
        return (2*x - (1))*(2*x - (1))
    elif (deg == 6):
        return (2*x - (1))*((8*x**2) - (8*x) + 1)
    elif (deg == 7):
        return (2*x - (1))*((32*x**3) - (48*x**2) + (18*x) - 1)
    elif (deg == 8):
        return ((8*x**2) - (8*x) + 1)
    elif (deg == 9):
        return ((8*x**2) - (8*x) + 1)*(2*x - (1))
    elif (deg == 10):
        return ((8*x**2) - (8*x) + 1)*((8*x**2) - (8*x) + 1)
    elif (deg == 11):
        return ((8*x**2) - (8*x) + 1)*((32*x**3) - (48*x**2) + (18*x) - 1)
    elif (deg == 12):
        return ((32*x**3) - (48*x**2) + (18*x) - 1)
    elif (deg == 13):
        return ((32*x**3) - (48*x**2) + (18*x) - 1)*(2*x - (1))
    elif (deg == 14):
        return ((32*x**3) - (48*x**2) + (18*x) - 1)*((8*x**2) - (8*x) + 1)
    elif (deg == 15):
        return ((32*x**3) - (48*x**2) + (18*x) - 1)*((32*x**3) - (48*x**2) + (18*x) - 1)
    else:
        return 0








def sumbasis1(X, n):
    n = n + 1
    suMatrix = []
    j = 0
    while j<n*n:
        sums = 0
        i = 0
        while i< len(X):
            sums = sums + (basis(j, X[i]))
            i = i + 1
        suMatrix.append(sums)
        j = j+1
    return suMatrix

def makematbasis(suMatrix, n):
    n = n + 1
    m = [[0 for i in range(n)]for j in range(n)]
    i = 0
    while i<n:
        j = 0
        while j<n:
            m[i][j] = suMatrix[(4*i)+(j%4)]
            j = j+1
        i = i + 1
    return m











def sumbasis2(X, Y, n):
    n = n+1
    suMatrix = []
    j = 0
    while j<n:
        sums = 0
        i = 0
        while i< len(X):
            sums = sums + (basis(j, X[i]))*Y[i]
            i = i + 1
        suMatrix.append(sums)
        j = j+1
    return suMatrix


def fitbasis(X,Y,degree):
    k = sumbasis1(X, degree)        
    m = makematbasis(k, degree)      
    qw = sumbasis2(X, Y, degree)
    X = gauss_jordan(m, qw)
    return X[0], X[1], X[2], X[3]



def cond_number(X, degree, sum=sum1, mat=makemat):
    k = sum(X, degree)
    m = mat(k, degree)
    return np.linalg.cond(m)
