import matplotlib.pyplot as plt
import library

def f(x, r=1):
    return 4*(r**2 - x**2)   #integral for volume of Steinmetz

print("For (i)")
print("Step \t \t Volume")
X = []
Y = []
for i in range(100,50100,100):
    pi = library.MonteCarlos(f,-1,1,i, 7865, 57657, 100)
    print(i ,"\t \t" ,pi)
    X.append(i)
    Y.append(pi)

plt.scatter(X,Y, s=1, color='red', label='Integration value at each step')
plt.axhline(y=Y[-1], label='Converges to '+str(round(Y[-1], 5)))
plt.xlabel("No. of steps")
plt.ylabel("Value of Integration")
plt.legend()
plt.grid()
plt.show()




















    
    
    


