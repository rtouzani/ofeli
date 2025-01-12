/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2025 Rachid Touzani
 
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

                      An example of a code using OFELI

   Solution of the 1-D linear transport equation using the Lax-Wendroff scheme

  ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   const double L = 10., c = 4.;

   banner();
   cout << "A PROGRAM TO TEST A 1-D LINEAR TRANSPORT EQUATION" << endl << endl;
   if (argc<3) {
      cout << "Usage: " << argv[0] << " <nx> <dt>" << endl;
      return EXIT_FAILURE;
   }
   ofstream pf("output.dat");

// Declare problem data
   try {
      size_t nx = atoi(argv[1]);
      double h = L/double(nx);
      theTimeStep = atof(argv[2]);
      theFinalTime = 1.;
      double cfl = c*theTimeStep/h;
      Vect<double> x, u(nx+1), v(nx+1);
      x.setUniform(0.,L,nx+1);

//    Set initial condition
      u.set("exp(-20*(x-2)^2)",x);

//    Time loop
      TimeLoop {
         v(1) = u(1) - cfl*(0.5*u(2)-cfl*(u(2)-2*u(1)));
         for (size_t i=2; i<=nx+1; i++)
            v(i) = u(i) - cfl*(0.5*(u(i+1)-u(i-1))-cfl*(u(i+1)-2*u(i)+u(i-1)));
         u = v;

//       save result in file each time step
         for (size_t i=1; i<=nx+1; i++)
            pf << x(i) << "   " << u(i) << endl;
         pf << endl;
      }

      cout << "Number of grid intervals: " << nx << endl;
      cout << "Time step:                " << theTimeStep << endl;
      cout << "CFL:                      " << cfl << endl;
   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}
