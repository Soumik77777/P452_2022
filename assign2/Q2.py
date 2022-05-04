import matplotlib.pyplot as plt
import library

def f(x):
    return 4*((1-x**2)**(1/2))

print("For (i)")
print("Step \t \t Pi")
X = []
Y = []
for i in range(1000,100000,100):
    pi = library.MonteCarlos(f,0,1,i, 65, 1021, 61)
    print(i ,"\t \t" ,pi)
    X.append(i)
    Y.append(pi)

plt.plot(X,Y,label="case 1")


print("For (ii)")
print("Step \t \t Pi")
X = []
Y = []
for i in range(1000,120000,100):
    pi = library.MonteCarlos(f,0,1,i, 572, 16381, 100)
    print(i ,"\t \t" ,pi)
    X.append(i)
    Y.append(pi)

plt.plot(X,Y,label="case 2")
plt.xlabel("No of Iteration")
plt.ylabel("Value of Integration")
plt.legend()
plt.grid()
plt.show()


print("In first choice of a and m, the value of pi oscillates around the real value and converges with a much higher frequency of oscilation.")
print("In second case, the oscillation happens with much lower frequency.")























    
    
    


