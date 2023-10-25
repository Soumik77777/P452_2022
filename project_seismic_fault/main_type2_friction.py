import math
import csv
import random
import numpy as np
import scipy
import scipy.special
import matplotlib.pyplot as plt

import fault_library

#######



def Func(sigma_by_P, tau_inf, friction, var):                   # var = [lambda, phi_j]
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
        func_vec[i,0] = - slip_term[i,0] + sigma_by_P - (tau_inf/friction) - fault_library.lambda_func(var[0], r[i])

    return func_vec


def solve_lam_phi(sigma_over_P, tau_inf, fric_p, init_guess, deg, tol=10e-6):
    r = fault_library.chebyshev_U_complimentary(deg)  ## coordinate
    s, w = fault_library.chebyshev_U_sw(deg)  ## s and w at each coordinate

    diff_sum = 10
    jk = 1  ## just to check no of iter
    while diff_sum > tol:
        jacobian_mat = np.zeros((deg + 1, deg + 1))
        for x in range(len(r)):
            jacobian_mat[x, 0] = fault_library.jaco_lambda(init_guess[0], r[x])
            for si in range(len(s)):
                jacobian_mat[x, si + 1] = - w[si] * fault_library.coeff_slip(r[x], s[si]) / (2 * np.pi)

        jacobian_inverse = np.linalg.inv(jacobian_mat)

        func_vec = Func(sigma_over_P, tau_inf, fric_p, init_guess)

        multiplied_term = np.matmul(jacobian_inverse, func_vec)

        new_var = []
        for i in range(len(init_guess)):
            new_var.append(init_guess[i] - multiplied_term[i, 0])

        diff_sum = abs((abs(new_var[0]) - abs(init_guess[0])) / init_guess[0])
        for i in range(len(new_var) - 1):
            diff_sum += abs((new_var[i + 1] - init_guess[i + 1]) / init_guess[i + 1])

        init_guess = new_var
        jk += 1

        print("Iteration ", str(jk))
        print("New variables=", new_var)
        print("diff= ", diff_sum)
        print()

    return init_guess


def delta_r(solution, deg, fric_p):              ## delta_bar multiplied by f
    soln_phi = solution[1:]
    del_r = []
    for ri in range(deg+1):
        del_r.append(0)
        for si2 in range(deg):
            if s[si2] > r[ri]:
                del_r[-1] += abs(soln_phi[si2] * w[si2]) * fric_p

    return del_r



def friction(f_p, f_r, delta_r, D_c):
    if delta_r>D_c:
        return f_r
    else:
        return f_p + ((f_r - f_p) * delta_r / D_c)


deg = 20
r = fault_library.chebyshev_U_complimentary(deg)  ## coordinate
s, w = fault_library.chebyshev_U_sw(deg)  ## s and w at each coordinate

lambda_0= 0.2
init_guess = [lambda_0]
for i in range(deg):
    init_guess.append(random.random())

sigma_over_P = 3.2
tau_inf = 3.2
fric_p = [1, 3, 5, 10, 15]
fric_r = fric_p[0]
D_c = 40           ## val of delta (x-axis scaling from -1 to 1)


solution = solve_lam_phi(sigma_over_P, tau_inf, fric_p[0], init_guess, deg=deg)
del_r = delta_r(solution, deg=deg, fric_p=fric_p[0])
plt.plot(r, del_r, '--', label='constant friction ($f_r$)= '+str(fric_p[0]))

solution = solve_lam_phi(sigma_over_P, tau_inf, fric_p[-1], init_guess, deg=deg)
del_r = delta_r(solution, deg=deg, fric_p=fric_p[-1])
plt.plot(r, del_r, '--', label='constant friction= '+str(fric_p[-1]))


all_delta =[]
for i in range(len(fric_p)):
    solution = solve_lam_phi(sigma_over_P, tau_inf, fric_p[i], init_guess, deg=deg)
    del_r = delta_r(solution, deg=deg, fric_p=fric_p[i])
    del_r_new = []
    for j in range(len(del_r)):
        if del_r[j]>D_c:
            soln2 = solve_lam_phi(sigma_over_P, tau_inf, fric_r, init_guess=init_guess, deg=deg)
            del_r_2 = delta_r(soln2, deg=deg, fric_p=fric_r)
            del_r_new.append(del_r_2[j])
        else:
            fr = friction(fric_p[i], fric_r, del_r[j], D_c)
            soln2 = solve_lam_phi(sigma_over_P, tau_inf, fr, init_guess=init_guess, deg=deg)
            del_r_2 = delta_r(soln2, deg=deg, fric_p=fr)
            del_r_new.append(del_r_2[j])
    all_delta.append(del_r_new)

    plt.plot(r, del_r_new, label='$f_p$='+str(fric_p[i]))

plt.xlabel("$r_i$", fontsize=14)
plt.ylabel("$\delta_r$", fontsize=14)

plt.legend()
plt.grid()
plt.show()







#### To plot r vs delta r
'''
soln_phi = solution[1:]
del_r = []
for ri in range(len(r)):
    del_r.append(0)
    for si2 in range(len(s)):
        if s[si2]>r[ri]:
            del_r[-1] += abs(soln_phi[si2] * w[si2]) * fric_p[0]

plt.plot(r, del_r, label='T= {con}'.format(con=str(round(fric_p[0],4))))
plt.plot([], [], label="n= 20, tol=10e-6")
plt.xlabel("$r_i$")
plt.ylabel("$\delta_r$")

plt.legend()
plt.grid()
plt.show()'''


'''
deg = 20

r = fault_library.chebyshev_U_complimentary(deg)        ## coordinate
s, w = fault_library.chebyshev_U_sw(deg)                ## s and w at each coordinate


##### constant with lambda
sigma_over_P = 3.2
tau_inf = 3.2
fric_p = [1]





for ik in range(len(fric_p)):
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



