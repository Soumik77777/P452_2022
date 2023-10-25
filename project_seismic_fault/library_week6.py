''' 
Personal library for useful functions; Soumik Bhattacharyya, Roll-1811155

This includes all necessary funcions used in P342,
but the ones that are required in the current week, are placed at the top.
'''

#import modules
import math
import csv
import random
#---------------


#Week6_Numerical Integration----------------------------------------------------

#Midpoint-----------------------------------------
def midpoint(a,b,N,f):
    h=(b-a)/N
    x=[]
    sum=0.0

    for i in range(N):
        xi=((a+i*h)+(a+((i+1)*h)))/2
        x.append(xi)

    for i in range(N):
        sum= sum + (h*f(x[i]))
    return sum

def midpoint_N(a,b,max_error,max_d2f): #This is for finding N for a given maximum error
    N = float(((((b-a)**3)*max_d2f)/(24*max_error))**(1/2))
    return math.ceil(N)


#Trapezoidal---------------------------------------
def trapezoidal(a,b,N,f):
    h = (b - a) / N
    x = []
    sum = 0.0

    for i in range(N+1):
        xi = a + (i * h)
        x.append(xi)

    for i in range(N+1):
        if i==0 or i==N:
            w=1
        else:
            w=2
        sum = sum + ((h*w*f(x[i]))/2)
    return(sum)

def trap_N(a,b,max_error,max_d2f): #This is for finding N for a given maximum error
    N = float(((((b-a)**3)*max_d2f)/(12*max_error))**(1/2))
    return math.ceil(N)


#simpson-------------------------------------------
def simpson(a,b,N,f):
    h = (b - a) / N
    x = [0 for i in range (N+1)]
    sum = 0.0

    for i in range(0,N+1,2):
        x[i] = a + (i*h)
    for i in range(1,N,2):
        x[i]=(x[i-1]+x[i+1])/2

    for i in range(N + 1):
        if i == 0 or i == N:
            w = 1
        elif i%2==0:
            w = 2
        else:
            w=4
        sum = sum + ((h*w*f(x[i]))/3)
    return sum

def simp_N(a,b,max_error,max_d4f): #This is for finding N for a given maximum error
    N = float(((((b-a)**5)*max_d4f)/(180*max_error))**(1/4))
    if N%2.0 == 0.0:
        return math.ceil(N)
    else:
        return math.ceil(N+1)


#Monte-carlo---------------------------------------
def monte_carlo(a,b,N,f):
    x=[]
    sum=0.0

    for i in range(N):
        r=random.uniform(a,b)    #random number in range (a,b)
        x.append(r)

    for i in range(N):
        sum+=f(x[i])

    value=((b-a)*sum)/N

    return value


#Week6 ends here-------------------------------------------------------------








#Functions for week5---------------------------------------------------------

def fix_ab (f, a, b):
    for i in range (1,10):
        if f(a)*f(b)<0:
            continue
        else:
            if abs(f(a))<abs(f(b)):
                a = a - 0.5*(b-a)
            else:
                b = b + 0.5*(b-a)
    return (a, b)

def bisection(f, a, b):
    a, b = fix_ab(f, a, b)

    print("Modified lower limit=", a)
    print("Modified upper limit=", b)

    e = 10 ** (-6)

    list_c = []

    while (b - a) > e:
        c = (a + b) / 2
        if f(c) == 0:
            print("The exact solution is, x=", c)
        else:
            if f(a) * f(c) < 0:
                a = a
                b = c
            else:
                a = c
                b = b
        list_c.append(c)

    print("The solution of the equation is, x=", c)

    return (c, list_c)

def list_error (list_c):
    error =[]
    for i in range(len(list_c)):
        error.append(abs(list_c[i]- list_c[i-1]))
        error[0] = "--"

    return (error)

def regula_falsi (f, a, b):
    a, b = fix_ab (f, a, b)

    print("Modified lower limit=", a)
    print("Modified upper limit=", b)

    e = 10 ** (-6)
    list_c = [a,b]
    while abs(list_c[-1]-list_c[-2])>e:
        c = (a*f(b) - b*f(a))/(f(b)- f(a))
        if f(c) == 0:
            print("The exact solution is, x=", c)
        else:
            if f(a) * f(c) < 0:
                a = a
                b = c
                list_c.append(round(c, 7))
            else:
                a = c
                b = b
                list_c.append(round(c, 7))

    print("The solution of the equation is, x=", c)

    del list_c[0]
    del list_c[1]

    return (c, list_c)


#Functions for Newton-Raphson

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

#functions for laguerre method

def poly_function(x,a):
    n=len(a)
    sum=0.0
    for i in range(n-1,-1,-1):
        sum+=a[i]*(x**i)
    return sum

def first_derivative(x,a):
    h=10**(-6)
    y=(poly_function(x+h,a)-poly_function(x-h,a))/(2*h)
    return y

def s_derivative(x,a):
    h = 10**(-6)
    y = (poly_function(x + h, a) + poly_function(x - h, a)-2*poly_function(x,a)) / (2 * h*h)
    return y

def deflate(sol, a):
    n = len(a)
    q = [0 for i in range(n - 1)]
    q[n - 2] = a[n - 1]
    for i in range(n - 3, -1, -1):
        q[i] = a[i + 1] + (sol * q[i + 1])

    return q


def laguerre(poly_function, first_derivative, s_derivative, a, i):
    n = len(a)
    if n != 2:
        j = i
        j1, j2 = i, 0
        k = 1
        if poly_function(i, a) != 0:
            while abs(j2 - j1) > 1E-4 and k < 200:
                g = first_derivative(j, a) / poly_function(j, a)

                h = g ** 2 - (s_derivative(j, a) / poly_function(j, a))
                f = ((n - 1) * (n * h - g ** 2))

                denominator1 = g + math.sqrt(f)

                denominator2 = g - math.sqrt(f)

                if abs(denominator1) > abs(denominator2):
                    j = n / denominator1
                else:
                    j = n / denominator2

                if k % 2 == 0:
                    j1 = j2 - j
                    j = j1

                else:
                    j2 = j1 - j
                    j = j2

                k += 1
        if k % 2 == 0:
            print(j1)
            a = deflate(j1, a)
        else:
            print(j2)
            a = deflate(j2, a)
        return a
    else:
        if a[1] * a[0] < 0 or a[1] < 0:
            print(a[0])
        else:
            print(-a[0])

        return 0



#Week5 ends here--------------------------------------------------------------------

def read_matrix (txt):
    with open(txt, 'r') as a:
        matrix = [[int(num) for num in row.split(' ')] for row in a]

    return matrix


def matrix_multiplication (matrix_m, matrix_n, n):
    matrix_L = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                matrix_L[i][j] += matrix_n[k][j] * matrix_m[i][k]
    return matrix_L


def gauss_jordan (m, v, n):
    partial_pivot(m, v, n)
    for i in range (n):
        factor1 = m[i][i]
        for j in range (n):
            m[i][j] = m[i][j]/factor1
        for q in range (n):
            v[i][q] = v[i][q]/factor1

        for k in range (n):
            if k!=i and m[k][i]!=0:
                factor2 = m[k][i]
                for l in range (i,n):
                    m[k][l] = m[k][l] - factor2*m[i][l]
                for r in range (n):
                    v[k][r] = v[k][r] - factor2* v[i][r]
    return (m,v)


def partial_pivot (m, v, n):
    for i in range (n-1):
        if m[i][i] ==0:
            for j in range (i+1,n):
                if abs(m[j][i]) > abs(m[i][i]):
                    m[i], m[j] = m[j], m[i]
                    v[i], v[j] = v[j], v[i]
    return (m,v)


def lu_decomposition (matrix, n):

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


def inverse_by_lu_decomposition (matrix, n):

    identity = [[0 for i in range(4)] for j in range(4)]
    for i in range(4):
        identity[i][i] = 1
    x = []

    '''The following process could be done easily using a for loop.
    But for those matrices, which are changed by partial pivoting, the columns of final inverse are interchanged.
    But pivoting is important, because at on estep in decomposition, it is divided by diagonal element
    So it is done manually for each row, and result of substitution is appended at the end of each step.'''

    matrix_0 = matrix.copy()
    partial_pivot(matrix_0, identity[0], n)
    (lower_0, upper_0) = lu_decomposition(matrix_0, n)
    x0 = forward_backward_substitution(lower_0, upper_0, identity[0], n)
    x.append(x0)

    matrix_1 = matrix.copy()
    partial_pivot(matrix_1, identity[1], n)
    (lower_1, upper_1) = lu_decomposition(matrix_1, n)
    x1 = forward_backward_substitution(lower_1, upper_1, identity[1], n)
    x.append(x1)

    matrix_2 = matrix.copy()
    partial_pivot(matrix_2, identity[2], n)
    (lower_2, upper_2) = lu_decomposition(matrix_2, n)
    x2 = forward_backward_substitution(lower_2, upper_2, identity[2], n)
    x.append(x2)

    matrix_3 = matrix.copy()
    partial_pivot(matrix_3, identity[3], n)
    (lower_3, upper_3) = lu_decomposition(matrix_3, n)
    x3 = forward_backward_substitution(lower_3, upper_3, identity[3], n)
    x.append(x3)

    inverse = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            inverse[i][j] = round(x[j][i], 3)

    return (inverse)

