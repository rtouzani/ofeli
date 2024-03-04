This directory contains examples of script files to run rita
in order to solve some ordinary differential equations

example1.rita:
Solution of the ODE:

y'(t) = y(t) + exp(t)  t>0,
y(0) = 0.

The solution is given by y(t) = t*exp(t)
We use the Heun (2nd order) scheme and compute error.
The result is stored in a file for plotting purposes.

example2.rita
Solution of an ODE giving a chaotic solution by the RK4 method

example3.rita:
Solution of a system of Lorentz differential system involving attractors.
We use the RK4 (4-th order Runge-Kutta) scheme
Solution is stored in file as well as phase portraits

