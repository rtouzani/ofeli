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

          A simple code to test the differential equation solver

    This code tests the ODE solver with the sample differential system:

                      A1*y' + A0*y = f

    where
              /          \         /         \       /                 \
              |  0   -1  |         | 1    0  |       |        0        |
         A0 = |          |,   A1 = |         |,  b = |                 |
              |  1   -1  |         | 0    1  |       | 3*(t-1)*exp(-t) |
              \          /         \         /       \                 /

   The solution is given by:
      y1 = t*exp(-t), y2 = (1-t)*exp(-t)
      
   Argument for execution: Time step
   You can test the code with various schemes by modifying the line:
        ODESolver ode(RK4,theTimeStep,theFinalTime,2);
   Just replace RK4 by FORWARD_EULER, BACKWARD_EULER, CRANK_NICOLSON,
   HEUN, AB2, BDF2

  ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   theFinalTime = 1.;
   if (argc<2) {
      cout << "Usage: ode_demo3 <time step>" << endl;
      return 0;
   }
   theTimeStep = atof(argv[1]);

// Solution as a system of 2 first-order ODEs
   try {
      ODESolver ode(BDF2,theTimeStep,theFinalTime,2);

//    Set differential equation coefficients
      DMatrix<double> A0(2,2), A1(2,2);
      ode.setMatrices(A0,A1);
      Vect<double> y(2), f(2), ff(2);
      y(1) = 0; y(2) = 1;
      f(1) = 0; f(2) = -3;
      ode.setInitial(y);
      ode.setInitialRHS(f); // Used for a multistep scheme only
      ode.setRHS(f);
      ode.setRK4RHS(ff);    // Used only if the RK4 scheme is used

//    Loop on time steps to run the time integration scheme
      TimeLoop {
         A1 = 1;
         A0(1,1) = 0; A0(1,2) = -1; A0(2,1) = 1; A0(2,2) = -1;
         f(1) = 0; f(2) = 3*(theTime-1)*exp(-theTime);

//       Vector ff might be used in the case of a multistep scheme like RK4
         double t=theTime-0.5*theTimeStep;
         ff(1) = 0; ff(2) = 3*(t-1)*exp(-t);

//       Run one time step
         ode.runOneTimeStep();
         cout << y;
      }

//    Output differential equation information and error
      cout << ode << endl;
      cout << "Error: " << fabs(y(1)-theFinalTime*exp(-theFinalTime)) << ", "
           << fabs(y(2)-(1-theFinalTime)*exp(-theFinalTime)) << endl;
   } CATCH_EXCEPTION

   return 0;
}
