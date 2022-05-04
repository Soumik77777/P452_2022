import library
import matplotlib.pyplot as plt
import numpy as np


M = library.constructM("assign2fit.txt", 21, 2)
i = 0
X = []
Y = []
while i < 21:
    X.append(M[i][0])
    Y.append(M[i][1])
    i = i + 1
plt.plot(X,Y,label="Data points")


a,b,c,d = library.fit(X,Y, 3)
print("Normal basis is coefficients: ",a,b,c,d)

def func(x):
    F = a + b*x + c*x**2 + d*x**3
    return F
A,B = library.Plot_Analytical(func, 0,1,0.001)
plt.plot(A,B, label = "Normal basis fit")

a,b,c,d = library.fitbasis(X,Y, 3)
print("Chebyshev basis coefficients: ",a,b,c,d)

def func2(x):
    F = a + (b*(2*x - (1))) + (c*((8*x**2) - (8*x) + 1)) + (d*((32*x**3) - (48*x**2) + (18*x) - 1))
    return F
A,B = library.Plot_Analytical(func2, 0,1,0.001)
plt.plot(A,B, label = "Chebyshev basis")
a,b,c,d = library.fit(X,Y, 3)

plt.legend()
plt.show()




print("Condition Number for Normal Basis",
      library.cond_number(X, 3, sum=library.sum1, mat=library.makemat))
print("Condition Number for Chebyshev Basis",
      library.cond_number(X, 3, sum=library.sumbasis1, mat=library.makematbasis))


print("The condition number for normal fitting is much higher than that of Chebyshev Basis."
      "Thus the fitting with Chebyshev basis is better thanthe normal one.")

'''
Normal basis is coefficients:  0.5746586674196168 4.725861442141896 -11.12821777764321 7.66867762290942
Chebyshev basis coefficients:  1.1609694790335525 0.39351446798815237 0.046849832090106576 0.23964617571596986
Condition Number for Normal Basis 12104.948671031543
Condition Number for Chebyshev Basis 3.8561465786155806
'''





    
    
    


