
import numpy as np
import matplotlib.pyplot as plt


a = 1 # lattice parameter

N = [i for i in range(1, 90)]
M = [i for i in range(1, 90)]

def L_supercell(n, m):
    return np.sqrt(n**2 + m**2 + n*m)

def theta_twisted(n, m):
    return np.degrees(np.arccos((n**2 + 4*n*m + m**2)/(2*(n**2 + n*m + m**2))))

data = []

for n in N:
    for m in M:
        L = L_supercell(n, m)
        theta = theta_twisted(n, m)
        data.append((n, m, L, theta))
        
        if abs(theta - 1.1) < 5e-2:
            print(f"n: {n}, m: {m}, L: {L:.2f} Angstroms, Theta: {theta:.4f} degrees")

data = np.array(data)        

plt.plot(data[:, 3], data[:, 2], 'o')
plt.xlabel('Theta (degrees)')
plt.ylabel('L (angstroms)')
plt.title('Twisted Bilayer WSe2')
plt.show()