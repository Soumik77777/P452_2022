import numpy as np
import matplotlib.pyplot as plt


gridpoints = 50
phi = np.zeros((gridpoints, gridpoints))
for row in range(len(phi)):
	phi[row][0] = 1

def laplace(phi, epsilon):
    d = 1.0
    while (d > epsilon):
        phi_next = phi.copy()
        phi[1:-1, 1:-1] = (phi_next[1:-1, :-2]+phi_next[1:-1, 2:]+phi_next[:-2, 1:-1]+phi_next[2:, 1:-1])/4
        d = np.sqrt(np.sum((phi-phi_next)**2)/np.sum(phi_next**2))
    return phi
V_l = laplace(phi, 1e-12)

##plotting phi
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
x, y = np.linspace(0, 1, gridpoints), np.linspace(0, 1, gridpoints)
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, V_l)
ax.set_xlabel('x-axis')
ax.set_ylabel('y-axis')
ax.set_zlabel('phi')
plt.show()

