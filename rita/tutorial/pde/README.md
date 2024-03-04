This directory contains examples of script files to run rita
in order to solve some partial differential equations

example1.rita:
Solution of a 1-D elliptic problem by P1 finite element method
An analytical solution is tested and error is given.

example2.rita:
Solution of the 2-D Laplace equation by P1 finite element method.
An analytical solution is tested and error is given.
The result is stored in a file for plotting using Gmsh.

example3.rita:
Solution of the 2-D heat equation by P1 finite element method.
The mesh is created step by step
Time integration uses the implicit Euler scheme
The result is stored in a file for plotting using Gmsh.

example4.rita:
Solution of a 2-D linear elasticity problem using the Q1 finite element method.
We consider a cantilever beam example.
The mesh file is given already.
The result is stored in a file for plotting using Gmsh.

example5.rita:
Solution of time dependent incompressible Navier-Stokes equations in 2-D
Numerical solution uses a 2nd-order projection time integration scheme.
Space discretization uses a stabilized P1/P1 finite element method
The numerical test concerns classical flow over a step.

