"""
1D Hydro Solver

@author: Samuel Gagnon-Hartman
26 February 2020
"""
import numpy as np
import matplotlib.pyplot as plt

# number of steps
nsteps = 1000

# set up parameters
dt = 1
n = 100
dx = 2.0
c_s = 1.0
x = np.arange(n) * dx
# f1 = np.ones(n)
# f2 = np.zeros(n)
J1 = np.zeros(n)
J2 = np.zeros(n)

A = 1.0
f1 = 1.0 + A * np.exp(-(x - 100) ** 2 / 40)
f2 = A * np.exp(-(x - 100) ** 2 / 40)

# plotting
plt.ion()
fig = plt.figure()
plt1, = plt.plot(x, f1, '-b')
plt.xlim([0, 200])
plt.ylim([0.5, 2.0])
fig.canvas.draw()

i = 0
while i < nsteps:
    # first consider no source term
    u = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))
    for i in range(n - 1):
        if u[i] > 0.0:
            J1[i] = f1[i] * u[i]
            J2[i] = f2[i] * u[i]
        else:
            J1[i] = f1[i + 1] * u[i]
            J2[i] = f2[i + 1] * u[i]

    # solve continuity equation
    f1[1:-1] = f1[1:-1] - (dt / dx) * (J1[1:-1] - J1[:-2])
    # solve boundaries
    f1[0] = f1[0] - (dt / dx) * J1[0]
    f1[-1] = f1[-1] + (dt / dx) * J1[-2]
    # solve Euler equation
    f2[1:-1] = f2[1:-1] - (dt / dx) * (J2[1:-1] - J2[:-2])
    # solve boundaries
    f2[0] = f2[0] - (dt / dx) * J2[0]
    f2[-1] = f2[-1] + (dt / dx) * J2[-2]

    # add in source term
    f2[1:-1] = f2[1:-1] - c_s**2 * (f1[2:] - f1[:-2]) / (2.0 * dx)
    f2[0] = f2[0] - 0.5 * c_s**2 * (f1[1] - f1[0]) / dx
    f2[-1] = f2[-1] - 0.5 * c_s**2 * (f1[-1] - f1[-2]) / dx
    f2[-1] = f2[-1] - 0.5 * c_s**2 * (f1[-1] - f1[-2]) / dx

    # update the plot
    plt1.set_ydata(f1)
    fig.canvas.draw()
    plt.pause(0.001)

    i += 1