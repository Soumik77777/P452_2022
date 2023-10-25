import math
import csv
import random
import numpy as np
import scipy
import scipy.special
import matplotlib.pyplot as plt

import fault_library

#######

def Func_axisym(var_phi):

    deg = len(var_phi)
    r = fault_library.chebyshev_U_complimentary(deg)  ## coordinate
    s, w = fault_library.chebyshev_T_sw(deg)  ## s and w at each coordinate

    slip_term = np.zeros((deg-1, deg))      ## coeff of phi terms multiplied by phi
    for x in range(deg - 1):
        for si in range(deg):
            slip_term[x, si] = - w[si] * fault_library.coeff_slip(r[x], s[si]) * var_phi[si] / (2 * np.pi)

    func_vector = np.zeros((deg, 1))  ## for i= 1 to n-1 eq from cheb, nth eq from interpolation
    for i2 in range(deg-1):
        sum2 = 1
        for j2 in range(deg):
            sum2 += slip_term[i2, j2]
        func_vector[i2][0] = sum2

    ## nth term
    sum3 = 0
    for i3 in range(1, deg+1):
        numer = np.sin((np.pi * (2 * deg - 1) * (2 * i3 - 1)) / (4 * deg))
        denom = np.sin((np.pi * (2 * i3 - 1)) / (4 * deg))
        sum3 += numer * var_phi[deg - i3] / denom                       ## +1 is omitted cz i3 is starting from 1

    func_vector[-1][0] = sum3/deg

    return func_vector


deg = 2

r = fault_library.chebyshev_T_complimentary(deg)        ## coordinate   (n-1) no
s, w = fault_library.chebyshev_T_sw(deg)                ## s and w at each coordinate

init_phi_j = []
for i4 in range(deg):
    init_phi_j.append(-3)

diff_sum = 10
iter_no = 0

while diff_sum>=0.1:
    print(iter_no)

    jacobian_mat = np.zeros((deg, deg))
    for i5 in range(deg-1):
        for i6 in range(deg):
            jacobian_mat[i5, i6] = - w[i6] * fault_library.coeff_slip(r[i5], s[i6]) / (2 * np.pi)

    for i7 in range(deg):
        sum4 = 0
        for i3 in range(1, deg + 1):
            numer = np.sin((np.pi * (2 * deg - 1) * (2 * i3 - 1)) / (4 * deg))
            denom = np.sin((np.pi * (2 * i3 - 1)) / (4 * deg))
            sum4 += numer * init_phi_j[deg - i3] / denom  ## +1 is omitted cz i3 is starting from 1
        jacobian_mat[-1, i7] = sum4 / deg

    jacobian_inverse = np.linalg.inv(jacobian_mat)

    func_vec = Func_axisym(init_phi_j)

    multiplied_term = np.matmul(jacobian_inverse, func_vec)

    new_var = []
    for i in range(len(init_phi_j)):
        new_var.append(init_phi_j[i] - multiplied_term[i, 0])

    diff_sum = 0
    for i in range(len(new_var)):
        diff_sum += abs((new_var[i] - init_phi_j[i]) / init_phi_j[i])

    init_phi_j = new_var
    iter_no += 1

    print("New variables=", new_var)
    print("diff= ", diff_sum)
    print()


'''

##### constant with lambda
constant_term = [0.2 * no for no in range(1)]
lambda_val = []
all_var_sol = []

for ik in range(len(constant_term)):
    lambda_0 = 0.2
    init_guess = [lambda_0]
    for i in range(deg):
        init_guess.append(random.random())

    diff_sum = 10
    jk = 1  ## just to check no of iter
    while diff_sum>0.000001:
        jacobian_mat = np.zeros((deg + 1, deg + 1))
        for x in range(len(r)):
            jacobian_mat[x, 0] = fault_library.jaco_lambda(init_guess[0], r[x])
            for si in range(len(s)):
                jacobian_mat[x, si+1] = w[si] * fault_library.coeff_slip(r[x], s[si]) / (2 * np.pi)

        jacobian_inverse = np.linalg.inv(jacobian_mat)

        func_vec = Func(constant_term[ik], init_guess)

        multiplied_term = np.matmul(jacobian_inverse, func_vec)

        new_var = []
        for i in range(len(init_guess)):
            new_var.append(init_guess[i] - multiplied_term[i,0])

        diff_sum = abs((abs(new_var[0]) - abs(init_guess[0])) / init_guess[0])
        print(diff_sum)
        for i in range(len(new_var)-1):
            diff_sum += abs((new_var[i+1] - init_guess[i+1]) / init_guess[i+1])


        init_guess = new_var
        jk += 1

        print("Value,", str(ik+1), "Iteration ", str(jk))
        print("New variables=", new_var)
        print("diff= ",diff_sum)
        print()

    lambda_val.append(abs(new_var[0]))
    all_var_sol.append(new_var)
'''



### Plot between lambda/ T and summation of F
'''
fun_values_sum = []
for i in range(len(constant_term)):
    func_values = Func(constant_term[i], all_var_sol[i])
    sum = 0
    for j in range(len(func_values)):
        sum += func_values[j]
    fun_values_sum.append(float(sum/(deg+1)))

#print(fun_values_sum)
#print(lambda_val)

#plt.title("Summation of functions as a function of lambda")
#plt.plot(lambda_val, fun_values_sum, label="n= 20, tol= 10e-6")
#plt.xlabel("$\lambda$")
plt.title("Summation of functions as a function of T")
plt.plot(constant_term, fun_values_sum, label="n= 20, tol= 10e-6")
plt.xlabel("T")
plt.ylabel("$\sum _i F_i$")

plt.legend()
plt.grid()
plt.show()
'''



### To plot Lambda vs T
'''
plt.plot(lambda_val, constant_term, label="Solution for lambda for different T")
plt.plot([], [], label="n= 200, tol= 10e-6")
plt.xlabel("Lambda")
plt.ylabel("T")

plt.legend()
plt.grid()
plt.show()
'''



#### To plot r vs delta r
'''
for i in range(len(all_var_sol)):
    soln_phi = all_var_sol[i][1:]
    del_r = []
    for ri in range(len(r)):
        del_r.append(0)
        for si2 in range(len(s)):
            if s[si2]>r[ri]:
                del_r[-1] += abs(soln_phi[si2] * w[si2])

    plt.plot(r, del_r, label='T= {con}'.format(con=str(round(constant_term[i],2))))
plt.plot([], [], label="n= 20, tol=10e-6")
plt.xlabel("$r_i$")
plt.ylabel("$\delta_r$")

plt.legend()
plt.grid()
plt.show()

'''


### To plot r vs del r
'''
coeff_mat = np.zeros((deg + 1, deg))        ## coeff of phi terms
for x in range(deg+1):
    for si in range(deg):
        coeff_mat[x, si] = w[si] * fault_library.coeff_slip(r[x], s[si]) / (2 * np.pi)

for i in range(len(all_var_sol)):
    soln_phi = all_var_sol[i][1:]
    del_r, del_r_tau = [], []
    for ri in range(len(r)):
        del_r.append(0)
        for si2 in range(len(s)):
            if s[si2]>r[ri]:
                del_r[-1] += abs(soln_phi[si2] * w[si2])
    for j in range(len(del_r)):
        sum = 0
        for k in range(len(s)):
            sum += abs(coeff_mat[j, k] * soln_phi[k])
        del_r_tau.append(del_r[j]/sum)

    plt.plot(r, del_r_tau, label='T= {con}'.format(con=str(round(constant_term[i],2))))

plt.plot([ ], [ ], label="n= 20, tol=10e-6")
plt.xlabel("$r_i$")
plt.ylabel("$\delta_r * \mu / \Delta Tau$")

plt.legend()
plt.grid()
plt.show()
'''



