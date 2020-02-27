# PHYS 432 Problem Set 4

A collection of simple simulations for PHYS 432 (Physics of Fluids) at McGill University.

## Prerequisites
This code runs on Python 3.6 or more recent versions and requires the ```numpy``` and ```matplotlib``` packages.

## File Descriptions

A usage guide and discussion for each of the python scripts included in this repository.

### advection.py

This script solves a 1D fluid using the FTCS method and Lax-Friedrich method separately. The results of using each method are plotted side-by-side as the script runs.

### adv_diff.py

This script solves the advection-diffusion equation using the Lax-Friedrich method for diffusion and the implicit method for advection. Running the script displays the evolution of a 1D fluid over time using these conditions.

### hydro.py

This script solves for the motion of a 1D fluid using the donor cell advection scheme with reflective boundary conditions. The conservative form of the hydro equations are used.

When the amplitude of the perturbation is increased, the main perturbation decays in amplitude quickly and the energy goes into activating lower modes in the fluid. Shock appears in the trail of the main wave. The width of the shock is set by the kinetic energy in the wave.

### hydro_bonus.py

The script expands on the script from ```hydro.py``` to now also
track the conservation of energy in the limit of a strong perturbation under an adiabatic process.

