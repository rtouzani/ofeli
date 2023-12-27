/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani
 
   This file is part of OFELI.

   OFELI is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   OFELI is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with OFELI. If not, see <http://www.gnu.org/licenses/>.

  ==============================================================================

                 An example of a Finite Element Code using OFELI

            Solution of a 1-D Parabolic problem using P1 Finite elements
            for space discretization and the Backward Euler scheme for 
            time discretization

  ==============================================================================*/

#include "OFELI.h"
#include "LinearPDE.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
// Read and output mesh data
   banner();
   cout << "A PROGRAM TO ILLUSTRATE A 1-D PARABOLIC EQUATION" << endl << endl;

// Read and set problem data
   theFinalTime = 1.0;
   if (argc<2) {
      cout << "Usage: " << argv[0] << " <nb of elements> [time step]" << endl;
      return EXIT_FAILURE;
   }
   int N = atoi(argv[1]);
   const double L = 1;
   Mesh ms(0.,L,N,true,1,1,1,1);
   theTimeStep = 0.1;
   if (argc>2)
      theTimeStep = atof(argv[2]);

   try {

//    Declare and initialize used vector
//    u: initial solution and solution at each time step
      Vect<double> u(ms);
      u.set("sin(2*pi*x)");
 
//    Instantiate equation class and declare used terms
      LinearPDE1D eq(ms);

//    Build the differential system
      TimeStepping ts(BACKWARD_EULER,theTimeStep,theFinalTime);
      ts.setPDE(eq);
      ts.setInitial(u);

//    Set solver of the linear system (See class LinearSolver for other choices)
      ts.setLinearSolver(CG_SOLVER,DILU_PREC);

//    Time loop
      TimeLoop {

//       Set right-hand side of equation and boundary conditions
         ts.setRHS("sin(2*pi*x)*exp(-t)*(4*pi*pi-1)");
         ts.setBC(1,"0.");

//       Set pde coefficients (=1) for the heat equation terms
         eq.setPDECoef(C10);
         eq.setPDECoef(C02);

//       Run the time step: The solution is stored in vector u
         ts.runOneTimeStep();
      }

//    Define exact solution and compute errors
      Vect<double> sol(ms);
      sol.setTime(theFinalTime);
      sol.set("sin(2*pi*x)*exp(-t)");
      cout << "Max-Norm error: " << (u-sol).Norm(NORM_MAX) << endl;

   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}