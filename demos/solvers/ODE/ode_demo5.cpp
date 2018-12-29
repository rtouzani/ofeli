/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

  Copyright (C) 1998 - 2019 Rachid Touzani

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

            Solution of the harmonic oscillator differential equation
                         using the ODESolver class

       We consider the ODE equation:
                     y'' + 0.15*y' + y = 0
       with the initial conditions: y(0) = 1, y(2) = 0

       We solve the problem by writing it as a system of 2 first-order ODE
       equations

  ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc<4) {
      cout << "Usage: ode_demo5 <time step> <final time> <output file>" << endl;
      return 0;
   }
   theTimeStep = atof(argv[1]);
   theFinalTime = atoi(argv[2]);
   ofstream outf(argv[3]);
   string ff = "phase.dat";
   if (argc>4)
      ff = argv[4];
   ofstream oph(ff.c_str());

// Solution is stored in file given by the fourth argument
// The phase portrait is solved in a file whose name can be given as a fifth argument,
// The default file name is 'phase.dat'

   try {
      ODESolver ode(FORWARD_EULER,theTimeStep,theFinalTime,2);

//    Set differential equation system
      ode.setF("y2");
      ode.setF("-0.15*y2-y1");
      Vect<double> y(2);
      y(1) = 1.; y(2) = 0.;
      ode.setInitial(y);

//    Loop on time steps to run the time integration scheme
      TimeLoop {
         ode.runOneTimeStep();
         outf << theTime << "  " << y(1) << endl;
         oph << y(1) << "  " << y(2) << endl;
      }

//    Output differential equation information and error
      cout << ode << endl;

   } CATCH_EXCEPTION

   return 0;
}
