/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

  Copyright (C) 1998 - 2020 Rachid Touzani

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

          A simple code to test the differential equation solver

    This code tests the ODE solver with the sample differential system:

                      [A2]{u''} + [A1]{u'} + [A0]{u} = {f}

    where
              /          \         /         \        /          \
              |  2   -1  |         | -2   0  |        |  2   0   |
         A0 = |          |,   A1 = |         |,  A2 = |          |,
              | -1    2  |         |  0   1  |        |  0   1   |
              \          /         \         /        \          /


              /                                                       \
             |   (2-2*pi^2)*sin(pi*t) - 2*pi*cos(pi*t) - cos(2*pi*t)  |
         f = |                                                        |.
             |  (2-4*pi^2)*cos(2*pi*t) - 2*pi*sin(2*pi*t) - sin(pi*t) |
              \                                                       /

   The solution is given by:
      u1 = sin(pi*t), u2 = cos(2*pi*t)
      
   Argument for execution: time step

  ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   theFinalTime = 1.;
   if (argc<2) {
      cout << "Usage: " << argv[0] << " <time step>" << endl;
      return 0;
   }
   theTimeStep = atof(argv[1]);
   Verbosity = 4;

// Solution as a system of second-order ODEs
   try {
      ODESolver ode(NEWMARK,theTimeStep,theFinalTime,2);

//    Set differential equation coefficients
      DMatrix<double> A0(2,2), A1(2,2), A2(2,2);
      Vect<double> u(2), v(2);
      u(1) = 0; u(2) = 1;
      v(1) = OFELI_PI; v(2) = 0;
      ode.setInitial(u,v);
      ode.setRHS("(2-2*pi^2)*sin(pi*t)-2*pi*cos(pi*t)-cos(2*pi*t)");
      ode.setRHS("(2-4*pi^2)*cos(2*pi*t)-2*pi*sin(2*pi*t)-sin(pi*t)");
      ode.setMatrices(A0,A1,A2);

//    Loop on time steps to run the time integration scheme
      TimeLoop {
//       We have to provide the matrices at each time step since in the case
//       of a direct solver, these matrices are modified after runOneTimeStep()
         A0(1,1) =  2.; A0(1,2) = -1.; A0(2,1) = -1.; A0(2,2) =  2.;
         A1(1,1) = -2.; A1(1,2) =  0.; A1(2,1) =  0.; A1(2,2) =  1.;
         A2(1,1) =  2.; A2(1,2) =  0.; A2(2,1) =  0.; A2(2,2) =  1.;
         ode.runOneTimeStep();
      }

//    Output differential equation information and error
      cout << ode << endl;
      cout << "Error: " << fabs(u(1)-sin(OFELI_PI*theFinalTime)) << ", "
                        << fabs(u(2)-cos(2*OFELI_PI*theFinalTime)) << endl;
   } CATCH_EXCEPTION

   return 0;
}
