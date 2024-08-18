import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 10.0  
q = 1000 
k = 200.0    
n = 10  # number of nodes

# Discretization
dx = L / (n - 1)
x = np.linspace(0, L, n)

# Matrix A and vector B for Fixed Temperature at Both Ends
T_0 = 0.0
T_L = 100.0
A = np.zeros((n, n))
B = np.zeros(n)



# Interior nodes
for i in range(1, n-1):
    A[i, i-1] = 1
    A[i, i] = -2 
    A[i, i+1] = 1
    B[i] = - q * dx**2 / k

# Boundary conditions
A[0, 0] = 1
B[0] = T_0
A[-1, -1] = 1
B[-1] = T_L
# Solve the linear system
T_fdm = np.linalg.solve(A, B)

# Analytical solution for Fixed Temperature at Both Ends
T_analytical = (T_L - T_0) / L * x + T_0 - (q / (2 * k)) * x**2 + (q * L*x) / (2 * k)

# Plotting results
plt.plot(x, T_fdm, 'o-', label='FDM Solution')
plt.plot(x, T_analytical, label='Analytical Solution')
plt.xlabel('Position along the rod (x)')
plt.ylabel('Temperature (T)')
plt.legend()
plt.title('1D Steady-State Heat Conduction with Fixed Temperature Boundary Conditions')
plt.show()
