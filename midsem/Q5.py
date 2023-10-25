###########

## Soumik Bhattacharyya, Roll- 1811155
## Q5. power iteration

import math

def modvec(x):
    suma = 0.0

    for i in range(len(x)):
        suma += (x[i]) * (x[i])
    return math.sqrt(suma)


def powermet(a, co, e=10e-5, ll=[], v1=[]):  # a is matrix, co is no of eigenvalues requested
    x = []  # ll and v1 are for storing eigenvalues and eigenvectors
    n = len(a)
    print("Please input guess matrix:")
    for i in range(n):
        x.append(float(input("i-th element:")))
    t = 1.0
    l = 1.0
    x1 = [0.0 for i in range(n)]
    if co > n:
        print("Function returned as no. of eigen value is greater than length of matrix")
        return

    while (t > e):

        # multiply, new x = A * old x
        for i in range(n):
            c = 0
            for j in range(n):
                c += a[i][j] * x[j]
            x1[i] = c

        # copy matrices new x to old x
        for i in range(n):
            x[i] = x1[i]

        # finding the max in x
        l1 = abs(x[0])
        for i in range(1, n):
            if abs(x[i]) > l1:
                l1 = abs(x[i])

        # dividing x by max element
        for i in range(n):
            x[i] = x[i] / l1

        # difference between old and new lambda
        t = abs(l1 - l)

        if t > e:
            l = l1  # if diff is more than tolerence, then old lambda=new lambda
            # the iteration continues
        else:
            print('Eigenvalue and corresponding eigenvector:')
            print(l1, x)  # print eigen value and corresponding eigenvector
            print()
            # storing eigenvalues and eigenvectors
            ll.append(l1)
            v1.append(x)
            co -= 1  # decrease eigenvalue no. counter
            u = []  # calculating normalised eigenvector U
            for i in range(n):
                u.append(x[i] / modvec(x))

            for i in range(n):  # calculating A* = A - lamda1 * U1 * transp(U)
                for j in range(n):
                    a[i][j] = a[i][j] - (l1 * u[i] * u[j])

            if co >= 1:  # continue finding next eigenvalue, if more than 1 is requested
                print("Next Eigenvalue---")
                powermet(a, co, e, ll, v1)
            return ll, v1


a,c=-1,-1
b = 2
n = 5

A=  [[2,  -1,   0,   0,   0],
    [-1,   2,  -1,   0,   0],
    [0,  -1,   2,  -1,   0],
    [0,   0,  -1,   2,  -1],
    [0,   0,   0,  -1,   2]]
print(A)
e=10e-5        # tolerence

print("No. of eigenvalues: (1/2)")
co=int(input())    # no of eigenvalues requested
lamb,vec= powermet(A,co,e)

print("Eigenvalues and corresponding eigenvectors:")
print(lamb,vec)

# checking if the first largest two eigenvalues and eigenvectors come from the given expression
# Comparison
print("Now, we compare the given expressions")
for k in range(1,3):   # k =1,2
    col,cov=0,0
    for i in range(1,len(vec[0])+1):
        # Calculating eigenvalues and eigenvector components from expressions
        lmd=b+2*math.sqrt(a*c)*math.cos((k*math.pi)/(n+1))
        vcc=2*((math.sqrt(c/a))**k)*math.sin((i*k*math.pi)/(n+1))
        print("Eigenvector component:",vcc)
        if round(vcc,2)==round(vec[k-1][i-1],2):
            print("verified.")
        else:
            print("Discrepency in eigen vector")
    print("Lambda 1:",lamb[k-1],", Lambda from the expression:",lmd)
    if round(lamb[k-1],2)==round(lmd,2): print("Lambda Verified")

