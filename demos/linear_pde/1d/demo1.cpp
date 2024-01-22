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

            Solution of a 1-D Elliptic problem using P1 Finite elements

  ==============================================================================*/

#include "OFELI.h"
#include "LinearPDE.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   const double L = 1;
   int N = 10;

// Read and output mesh data
   banner();
   cout << "A PROGRAM TO ILLUSTRATE A 1-D ELLIPTIC EQUATION" << endl << endl;
   if (argc>1)
      N = atoi(argv[1]);
   Mesh ms(0.,L,N,true,1,1,1,1);
   Vect<double> u(ms), f(ms), bc(ms);

   try {
//    Declare problem data (matrix, rhs, boundary conditions, body forces)
      f.set("(16*pi*pi+x)*sin(4*pi*x)");
      LinearPDE1D eq(ms,u);
      eq.setPDECoef(PDECoefType::C02,1.);
      eq.setPDECoef(PDECoefType::C00,"x");
      eq.setDirichlet(bc);
      eq.setBodyForce(f);

//    Solve problem
      eq.run();

//    Output solution and error
      cout << "\nSolution:\n" << u;
      Vect<double> sol(ms);
      sol.set("sin(4*pi*x)");
      cout << "Max-Norm error: " << (u-sol).Norm(NORM_MAX) << endl;
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}