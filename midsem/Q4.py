###########

## Soumik Bhattacharyya, Roll- 1811155
## Q4. The one on decaying fitting

import os
import csv
import math
import numpy as np
import matplotlib.pyplot as plt




###########

def read_csv(filename, location=None, poprows=None, delimiter=None):

    if location==0:
        filepath = filename
    elif location!= None and location!=0:
        filepath = str(location) + str(filename)
    else:
        filepath = str(os.getcwd()) + str('\\') + str(filename)

    if delimiter=='\t':
        delim = '\t'
    else:
        delim = ','

    with open(filepath, 'r') as infile:
        data = csv.reader(infile, delimiter=delim)
        datalist, rows = [], []
        for row in data:
            datalist.append(row)
        if poprows!=None:
            for i in range(poprows):
                datalist.pop(0)
        for j in range(len(datalist[0])):
            globals()['string%s' % j] = []
            for k in datalist:
                globals()['string%s' % j].append(float(k[j]))
            rows.append(globals()['string%s' % j])
        infile.close()

    return rows


data = read_csv('msfit.txt', delimiter='\t')
time_lst = data[0]
count_lst = data[1]
sigma_lst = data[2]


def chi2_linear_regression(xvals, yvals, sigma_lst):
    
    n = len(xvals)
    s, sx, sy, sxx, sxy = 0, 0, 0, 0, 0
    for i in range(n):
        s += 1 / sigma_lst[i]**2
        sx += xvals[i] / sigma_lst[i]**2
        sy += yvals[i] / sigma_lst[i]**2
        sxx += xvals[i]**2 / sigma_lst[i]**2
        sxy += xvals[i]*yvals[i] / sigma_lst[i]**2

    delta = s*sxx - sx**2
    a = (sxx*sy - sx*sxy) / delta
    b = (s*sxy - sx*sy) / delta

    dof = n - 2
    
    chi2 = 0
    for i in range(n):                          # chi-sqr calculation
        chi2 += (yvals[i] - a - b*xvals[i])**2 / sigma_lst[i]**2

    delA2 = sxx / delta; delB2 = s / delta
    cov = -sx / delta
    return a, b, delA2, delB2, cov, chi2


log_count_lst, new_sigma_lst = [], []
for i in range(len(count_lst)):
    log_count_lst.append(math.log(count_lst[i]))
    new_sigma_lst.append(1/sigma_lst[i])


params = chi2_linear_regression(xvals=time_lst, yvals=log_count_lst, sigma_lst=new_sigma_lst)

print(params)
print("Intercept: ", params[0])
print("Slope: ", params[1])
print("Chi square: ", params[-1])


calc_val = []
for i in time_lst:
    calc_val.append(i*params[1] + params[0])

plt.errorbar(time_lst, log_count_lst, yerr=new_sigma_lst, label='data points')
plt.plot(time_lst, calc_val, label='fitted straigh line')
plt.xlabel("Time Axis (sec)")
plt.ylabel("log(N)")

plt.legend()
plt.grid()
plt.show()


def half_life(a0, k):
    return a0 / (2 * abs(k))

'''print()
print("Half-life: ", half_life(params[0], params[1]), "seconds")
'''
print()
print("Average Life-time (1/rate constant): ", 1/abs(params[1]))