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
F1 = np.ones(n)
F2 = np.copy(F1)
F1[0] = 0.0
F2[0] = 0.0
D = 1.0
u = -0.01

# group useful constants
alpha = u * dt / dx
beta = D * dt / dx ** 2

# plotting
plt.ion()
fig, axes = plt.subplots(1, 2)
plt1, = axes[0].plot(x, F1, 'ro')
plt2, = axes[1].plot(x, F2, 'ro')
fig.canvas.draw()

count = 0
while count < nsteps:

    # constant boundaries
    F1[-1] = 1.0
    F1[0] = 0.0
    F2[-1] = 1.0
    F2[0] = 0.0

    # Lax-Friedrich advection
    F1[1:n - 1] = 0.5 * (F1[:n - 2] + F1[2:]) - 0.5 * alpha * (-F1[:n - 2] + F1[2:])

    # implicit diffusion
    A = np.eye(n) * (1 + 2 * beta) + np.eye(n, k=1) \
        * (-beta) + np.eye(n, k=-1) * (-beta)
    # enforce boundary conditions
    A[n - 1][n - 1] = 1.0
    A[n - 1][n - 2] = 0.0
    A[0][0] = 1.0
    A[0][1] = 0.0
    # solve for F2
    F2 = np.linalg.solve(A, F2)

    # update the plot
    plt1.set_ydata(F1)
    plt2.set_ydata(F2)
    fig.canvas.draw()
    plt.pause(0.001)
    count += 1
