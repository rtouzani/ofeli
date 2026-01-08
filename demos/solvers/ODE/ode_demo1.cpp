/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

  Copyright (C) 1998 - 2026 Rachid Touzani

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
                 Case of data given by algebraic expressions
 
   The tested differential equation is:
          y'(t) + y(t) = 2*exp(t)-1,   t>0
          y(0) = 0

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
      cout << "Usage: " << argv[0] << " <time step>" << endl;
      return EXIT_FAILURE;
   }
   theTimeStep = atof(argv[1]);

// Declare equation, give initial solution and give coefficients defining the ode
// Use the Forward Euler scheme (you can modify this)
   try {
      ODESolver ode(FORWARD_EULER);

//    Give initial condition
      ode.setInitial(0.);

//    Set function that defines the ODE
      ode.setF("2*exp(t)-1-y");

//    Solve the equation: Contains the loop on time steps
      ode.run();

//    Output differential equation information, numerical solution and error at
//    final time
      cout << ode << endl;
      cout << "Solution:  " << ode.get() << endl; 
      cout << "Error:     " << fabs(exp(theFinalTime)-1-ode.get()) << endl;
   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}