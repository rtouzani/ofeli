
This folder contains a list of simple demo programs in 1D using OFELI.
The library ofeli must have been installed before running these tests.

1. An elliptic 2-point value problem: elliptic

A demo program to solve a 1D elliptic equation and compute discretization error.
To execute the program, type

          elliptic <n>

where n is the number of finite elements. Its default value is 10

- A heat equation: heat

A demo program to solve the 1D heat equation using a 3-stencil finite difference scheme in space
and the Backward Euler scheme in time.
To execute the program, type

          heat <nx> <dt>

where n is the number of finite elements and dt is the time step.

- A linear transport equation: transport

A demo program to solve a 1D linear transport equation using a simple finite differnce upwind scheme
in space and the Backward Euler scheme in time
To execute the program, type

          transport <nx> <dt>

where nx is the number of space discretization points and dt is the time step.

As an example, here is the command line to create the executable
'elliptic':

         g++ -std=c++1y -I$CMAKE_INSTALL_PREFIX/include/ofeli elliptic.cpp -lofeli

where $CMAKE_INSTALL_PREFIX is the defined prefix path (or the default one).
