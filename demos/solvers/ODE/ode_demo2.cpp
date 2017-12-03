/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

  Copyright (C) 1998 - 2018 Rachid Touzani

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

                     A program to illustrate the ODE solver
                         Case of data given numerically
 
   The tested differential equation is:
          y'(t) + y(t) = 2*exp(t)-1,   t>0
          y(0) = 0.

   The exact solution is given by
          y(t) = exp(t) - 1

 ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
// Read or set some data (Final time and time step)
   theFinalTime = 1.;
   if (argc<2) {
      cout << "Usage: ode_demo2 <time step>" << endl;
      return 0;
   }
   theTimeStep = atof(argv[1]);

// Declare equation, give initial solution
// Use the Heun scheme (you can modify this)
   ODESolver ode(HEUN);
   ode.setInitial(0.);
   ode.setInitialRHS(1.);

/* Solve the equation
   The coefficients of the ode are given for each time step
   We also give the appropriate right-hand side for the RK4 scheme when this
   one is used                                                                 */
   TimeLoop {
      ode.setCoef(1,1,0,2*exp(theTime)-1);
      ode.setRK4RHS(2*exp(theTime-0.5*theTimeStep)-1);
      ode.runOneTimeStep();
   }

// Output differential equation information, numerical solution and error at
// final time
   cout << ode << endl;
   cout << "Error: " << fabs(exp(theFinalTime)-1-ode.get()) << endl;
   return 0;
}
