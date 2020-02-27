"""
Advection

@author: Samuel Gagnon-Hartman
25 February 2020
"""

import numpy as np
import matplotlib.pyplot as plt

# steps up the grid and flow parameters
n = 100
f_0 = np.linspace(0, 1, length)
f_1 = np.copy(f_0)
f_2 = np.copy(f_0)
x = np.linspace(0, 100, length)
nsteps = 2000
dx = 1
dt = 1
u = -0.1
step = u * dt / (2 * dx)

# set up the plot
plt.ion()
fig, axes = plt.subplots(1, 2)
axes[0].plot(x, f_0, 'k-')
axes[1].plot(x, f_2, 'k-')
plt1, = axes[0].plot(x, f_1, '.r')
plt2, = axes[1].plot(x, f_2, '.r')
fig.canvas.draw()

i = 0
while i < nsteps:
    # FTCS method
    f_1[1:n-1] = f_1[1:n-1] - step * (f_1[2:] - f_1[:n-2])
    # Lax-Friedrich method
    f_2[1: n-1] = 0.5*(f_2[2:] + f_2[:n-2]- step*(f_2[2:] - f_2[:n-2]))
    # plotting
    plt1.set_ydata(f_1)
    plt2.set_ydata(f_2)
    fig.canvas.draw()
    plt.pause(0.00001)
    i += 1