/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

  Copyright (C) 1998 - 2016 Rachid Touzani

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
 
                                t s _ d e m o _ 1
 
  A program to show the usage of the class TimeStepping for a time dependent problem

  The program solves the heat equation in 2-D using P1 finite elements in space
  and can use any of the implemented schemes for time integration.

  The choice of the diffusivity coefficient equal to 1 and the source term
       f(x,y,t) = -3*exp(-t)*exp(x+y)
  leads to the solution:
       u(x,y,t) = exp(-t)*exp(x+y)

  ==============================================================================*/

#include "OFELI.h"
#include "Therm.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
// Read and set problem data
   theFinalTime = 1.;
   if (argc<2) {
      cout << "ts_demo1 <mesh_file> [dt]" << endl;
      exit(1);
   }
   Mesh ms(argv[1],true);
   theTimeStep = 0.1;
   if (argc>2)
      theTimeStep = atof(argv[2]);

// Declare and initialize used vector
// u: initial solution and solution at each time step
   Vect<double> u(ms), f(ms), bc(ms);
   u.set("exp(x+y)");

// Instantiate equation class and declare used terms
   DC2DT3 eq(ms);
   eq.setTerms(CAPACITY|DIFFUSION|SOURCE);

// Build the differential system
// We use for instance the BDF2 scheme. Other choices are: FORWARD_EULER, BACKWARD_EULER,
// CRANK_NICOLSON, AB2 (Adams-Bashforth), HEUN, RK4, ...
   TimeStepping ts(BACKWARD_EULER,theTimeStep,theFinalTime);
   ts.setPDE(eq);
   ts.setInitial(u);
// Set solver of the linear system (See class LinearSolver for other choices)
   ts.setLinearSolver(CG_SOLVER,DILU_PREC);

// Time loop
   TimeLoop {

//    To give source term by an expression, set time value
      f.setTime(theTime);
      f.set("-3*exp(-t)*exp(x+y)");
      ts.setRHS(f);

//    Set Dirichlet boundary condition
      bc.setTime(theTime);
      bc.setNodeBC(1,"exp(-t)*exp(x+y)");
      ts.setBC(bc);

//    Run the time step: The solution is stored in vector u
      ts.runOneTimeStep();
   }

// Display TimeStepping class information
   cout << ts;

// Compute error
   Vect<double> U(ms,"u_ex",theFinalTime);
   U.set("exp(-1)*exp(x+y)");
   cout << "Solution L2-Norm: " << u.getWNorm2() << endl;
   cout << "Error in L2-Norm: " << (u-U).getWNorm2() << endl;
   return 0;
}
