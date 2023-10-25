import math
import csv
import random
import numpy as np
import scipy
import scipy.special
import matplotlib.pyplot as plt

import fault_library

#######



def Func(k, var):                   # k = constant term for friction, lam = length ratio
    deg = len(var) - 1

    func_vec = np.zeros((deg+1, 1))

    r = fault_library.chebyshev_U_complimentary(deg)  ## coordinate
    s, w = fault_library.chebyshev_U_sw(deg)  ## s and w at each coordinate

    coeff_mat = np.zeros((deg + 1, deg))        ## coeff of phi terms
    for x in range(deg+1):
        for si in range(deg):
            coeff_mat[x, si] = w[si] * fault_library.coeff_slip(r[x], s[si]) / (2 * np.pi)

    slip_term = np.zeros((deg + 1, 1))          ## phi multiplied with coeff
    for i in range(len(coeff_mat)):
        sum = 0
        for j in range(len(coeff_mat[0])):
            sum+= var[j+1]*coeff_mat[i, j]
        slip_term[i, 0] = sum


    for i in range(len(func_vec)):              ## final value of function
        func_vec[i,0] = - slip_term[i,0] + k - fault_library.lambda_func(var[0], r[i])

    return func_vec


deg = 20

r = fault_library.chebyshev_U_complimentary(deg)        ## coordinate
s, w = fault_library.chebyshev_U_sw(deg)                ## s and w at each coordinate


##### constant with lambda
constant_term = [0.5 * no for no in range(1, 7)]
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
                jacobian_mat[x, si+1] = - w[si] * fault_library.coeff_slip(r[x], s[si]) / (2 * np.pi)

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




### Plot between lambda/ T and summation of F

fun_values_sum = []
for i in range(len(constant_term)):
    func_values = Func(constant_term[i], all_var_sol[i])
    sum = 0
    for j in range(len(func_values)):
        sum += func_values[j]
    fun_values_sum.append(float(sum/(deg+1)))

#print(fun_values_sum)
#print(lambda_val)

plt.title("Summation of functions as a function of lambda")
plt.plot(lambda_val, fun_values_sum, label="n= 20, tol= 10e-6")
plt.xlabel("$\lambda$")
plt.ylabel("$\sum _i F_i / (n+1)$")

plt.legend()
plt.grid()
plt.show()

plt.title("Summation of functions as a function of T")
plt.plot(constant_term, fun_values_sum, label="n= 20, tol= 10e-6")
plt.xlabel("T")
plt.ylabel("$\sum _i F_i / (n+1)$")

plt.legend()
plt.grid()
plt.show()




### To plot Lambda vs T

plt.title("Solving Lambda for different T")
plt.plot(lambda_val, constant_term, label="n= 200, tol= 10e-6")
plt.xlabel("Lambda")
plt.ylabel("T")

plt.legend()
plt.grid()
plt.show()




#### To plot r vs delta r

for i in range(len(all_var_sol)):
    soln_phi = all_var_sol[i][1:]
    del_r = []
    for ri in range(len(r)):
        del_r.append(0)
        for si2 in range(len(s)):
            if s[si2]>r[ri]:
                del_r[-1] += abs(soln_phi[si2] * w[si2]) / constant_term[i]

    plt.plot(r, del_r, label='T= {con}'.format(con=str(round(constant_term[i],4))))
plt.xlabel("$r_i$")
plt.ylabel("$\delta_r$")
plt.title("")

plt.legend()
plt.grid()
plt.show()




### To plot r vs del r

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


plt.xlabel("$r_i$")
plt.ylabel("$\delta_r * \mu / \Delta Tau$")

plt.legend()
plt.grid()
plt.show()




