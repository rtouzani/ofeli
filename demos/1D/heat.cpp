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

       Solution of the 1-D heat equation using finite differences in space
                     and the Backward Euler scheme in time

  ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   const double L = 1.;

/// Read and output mesh data
   banner();
   cout << "A PROGRAM TO TEST THE HEAT EQUATION" << endl << endl;
   if (argc<3) {
      cout << "Usage: " << argv[0] << " <nx> <dt>" << endl;
      return EXIT_FAILURE;
   }
   ofstream pf("output.dat");

// Declare problem data (matrix, rhs, boundary conditions, body forces)
   try {
      size_t nx = atoi(argv[1]);
      double h = L/double(nx);
      theTimeStep = atof(argv[2]);
      theFinalTime = 1.;
      TrMatrix<double> A(nx+1);
      Vect<double> x, f(nx+1), u(nx+1), v(nx+1);
      x.setUniform(0.,L,nx+1);
      u.set("sin(pi*x)",x);

//    Time loop
      TimeLoop {

//       Build linear system
         f.setTime(theTime);
         f.set("exp(t)*sin(pi*x)*(1+pi^2)",x);
         for (size_t i=2; i<=nx; i++) {
            A(i,i  ) =  1./theTimeStep + 2./(h*h);
            A(i,i+1) = -1./(h*h);
            A(i,i-1) = -1./(h*h);
            v(i) = u(i)/theTimeStep + f(i);
         }

//       Impose Dirichlet boundary conditions at both ends
         A(1,1) = 1./theTimeStep;
         A(1,2) = 0.;
         v(1) = 0.;
         A(nx+1,nx+1) = 1./theTimeStep;
         A(nx+1,nx) = 0.;
         v(nx+1) = 0.;

//       Solve
         A.solve(v);
         u = v;

//       save result in file each time step
         for (size_t i=1; i<=nx+1; i++)
            pf << x(i) << "   " << u(i) << endl;
         pf << endl;
      }

      Vect<double> sol(nx+1);
      sol.setTime(theFinalTime);
      sol.set("exp(t)*sin(pi*x)",x);
      cout << "Number of grid points: " << nx << endl;
      cout << "Time step:             " << theTimeStep << endl;
      cout << "L2-Norm Error =        " << (u-sol).getWNorm2() << endl;
   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}
