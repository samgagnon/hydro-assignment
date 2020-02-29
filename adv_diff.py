"""
Advection-Diffusion Equation

@author: Samuel Gagnon-Hartman
26 February 2020
"""

import numpy as np
import matplotlib.pyplot as plt

nsteps = 1000  # number of steps

# initialize parameters
n = 100
dx = 1.0
dt = 1
x = np.linspace(0, n, 100)
F1 = np.linspace(0, 1, 100)
F2 = np.copy(F1)
F1[0] = 0.0
F2[0] = 0.0
D1 = 1.0
D2 = 10
u = -0.01

# group useful constants
alpha = u * dt / dx

beta1 = D1 * dt / dx ** 2
beta2 = D2 * dt / dx **2

# plotting
plt.ion()
fig, axes = plt.subplots(1, 2)
axes[0].set_title("D=1", fontsize=15)
axes[1].set_title("D=10", fontsize=15)
plt1, = axes[0].plot(x, F1, 'ro')
plt2, = axes[1].plot(x, F2, 'ro')
fig.canvas.draw()

i = 0
while i < nsteps:
    # D1    
    # implicit diffusion
    A1 = np.eye(n) * (1 + 2 * beta1) + (-beta1) * np.eye(n, k=1) + (-beta1) * np.eye(n, k=-1)
    # enforce boundary conditions
    A1[0][0] = 1.0
    A1[0][1] = 0.0
    A1[-1][-1] = 1.0 + beta1
    # Lax-Friedrich advection
    F1[1:n - 1] = 0.5 * (F1[:n - 2] + F1[2:]) - 0.5 * alpha * (-F1[:n - 2] + F1[2:])
    
    F1 = np.linalg.solve(A1, F1)
    
    # D2
    # implicit diffusion
    A2 = np.eye(n) * (1 + 2 * beta2) + (-beta2) * np.eye(n, k=1) + (-beta2) * np.eye(n, k=-1)
    # enforce boundary conditions
    A2[0][0] = 1.0
    A2[0][1] = 0.0
    A2[-1][-1] = 1.0 + beta2
    # solve for F2
    F2 = np.linalg.solve(A2, F2)
    
    F2[1:n - 1] = 0.5 * (F2[:n - 2] + F2[2:]) - 0.5 * alpha * (-F2[:n - 2] + F2[2:])
    # update the plot
    plt1.set_ydata(F1)
    plt2.set_ydata(F2)
    fig.canvas.draw()
    plt.pause(0.001)
    i += 1
