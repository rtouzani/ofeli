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

          Solution of the Lorenz dynamical system by the RK4 scheme
                       using the ODESolver class

       We consider the ODE system:

           y1' = 10*(y2-y1)
           y2' = y1*(27-y3) - y2
           y3' = y1*y2 - 8/3*y3

       with the initial condition y1 = 1, y2 = 0, y3 = 0

  ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc<4) {
      cout << "Usage: ode_demo7 <time step> <final time> <output file>" << endl;
      return 0;
   }
   theTimeStep = atof(argv[1]);
   theFinalTime = atoi(argv[2]);

// Solution is stored in file given by the fourth argument
// The phase portrait is solved in a file whose name can be given as a fifth argument,
// The default file name is 'phase.dat'

   ofstream outf(argv[3]);
   string ff = "phase.dat";
   if (argc>4)
      ff = argv[4];
   ofstream oph(ff);

// Solution as a system of 3 first-order ODEs
   try {
      ODESolver ode(RK4,theTimeStep,theFinalTime,3);

//    Set differential equation system
      ode.setF("10*(y2-y1)");
      ode.setF("y1*(28-y3)-y2");
      ode.setF("y1*y2-8./3.*y3");
      Vect<double> y(3), f(3), ff(3);
      y(1) = 1.; y(2) = 0., y(3) = 0.;
      ode.setInitial(y);

//    Loop on time steps to run the time integration scheme
      TimeLoop {

//       Run one time step
         ode.runOneTimeStep();
         outf << theTime << "  " << y(1) << "  " << y(2) << "  " << y(3) << endl;
         oph << y(1) << "  " << ode.getTimeDerivative(1) << "  " << ode.getTimeDerivative(2)
             << "  " << ode.getTimeDerivative(3) << "  " << endl;
      }

//    Output differential equation information and error
      cout << ode << endl;

   } CATCH_EXCEPTION

   return 0;
}
