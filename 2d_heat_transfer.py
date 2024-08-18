import numpy as np
import matplotlib.pyplot as plt

# Parameters
Lx = 1.0  
Ly = 1.0  
Nx = 80   
Ny = 80   
dx = Lx / (Nx - 1)  
dy = Ly / (Ny - 1)  
k = 15.0  
h = 10
T_inf = 30 
q = 1000  

# Create grid
T = np.zeros((Nx, Ny))

# Construct the coefficient matrix A and the right-hand side vector b
A = np.zeros((Nx * Ny, Nx * Ny))
b = np.zeros(Nx * Ny)

# Helper function to convert (i, j) index to linear index
def index(i, j):
    return i + j * Nx

for i in range(Nx):
    for j in range(Ny):
        idx = index(i, j)

        # Internal nodes
        if 0 < i < Nx - 1 and 0 < j < Ny - 1:
            A[idx, index(i-1, j)] = k / dx**2
            A[idx, index(i+1, j)] = k / dx**2
            A[idx, index(i, j-1)] = k / dy**2
            A[idx, index(i, j+1)] = k / dy**2
            A[idx, idx] = -2 * (k / dx**2 + k / dy**2)
            b[idx] = -q

        # Boundary conditions
        if i == 0:
            A[idx, idx] = 1
            if Nx > 1:
                A[idx, index(i+1, j)] = -k / (k+h*dx)
            b[idx] = (h*dx) * T_inf / (k + h * dx)

        elif i == Nx - 1:
            A[idx, idx] = 1
            if Nx > 1:
                A[idx, index(i-1, j)] = -k /(k+h*dx)
            b[idx] = (h* dx) * T_inf / (k + h * dx)

        elif j == 0:
            A[idx, idx] = 1
            if Ny > 1:
                A[idx, index(i, j+1)] = -k / (k+h*dy)
            b[idx] = (h* dy) * T_inf / (k + h * dy)

        elif j == Ny - 1:
            A[idx, idx] = 1
            if Ny > 1:
                A[idx, index(i, j-1)] = -k / (k+h*dy)
            b[idx] = (h* dy) * T_inf / (k + h * dy)

# Solve the linear system
T_flat = np.linalg.solve(A, b)

# Reshape the solution back to 2D
T = T_flat.reshape((Nx, Ny))

# Plot the result
plt.figure(figsize=(8, 6))
plt.contourf(T, cmap='hot', levels=50)
plt.colorbar(label='Temperature (°C)')
plt.title('Temperature Distribution in the 2D Block')
plt.xlabel('x direction (nodes)')
plt.ylabel('y direction (nodes)')
plt.show()

# Calculate total heat generated
total_heat_generated = q * dx * dy * Nx * Ny

# Calculate total heat lost through convection
heat_loss_x0 = np.sum(h * dx * (T[0, :] - T_inf))
heat_loss_xL = np.sum(h * dx * (T[-1, :] - T_inf))
heat_loss_y0 = np.sum(h * dy * (T[:, 0] - T_inf))
heat_loss_yL = np.sum(h * dy * (T[:, -1] - T_inf))

total_heat_lost = heat_loss_x0 + heat_loss_xL + heat_loss_y0 + heat_loss_yL
print('heat lost from simulation is',total_heat_lost)
print('heat generated is',q*Lx*Ly)