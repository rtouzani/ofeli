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

            Solution of a 2-D elliptic problem using P1 Finite elements

  ==============================================================================*/

#include "OFELI.h"
#include "LinearPDE.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
// Expand arguments
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <project_file>\n";
      return EXIT_FAILURE;
   }
   cout << "A PROGRAM TO ILLUSTRATE A 2-D ELLIPTIC EQUATION" << endl << endl;
   IPF data(string(argv[0])+" - 1.0",argv[1]); 

   try {

//    Define problem data
      Mesh ms(data.getMeshFile());
      Vect<double> bc(ms), bf(ms), u(ms);
      Prescription p(ms,data.getPrescriptionFile());
      p.get(BOUNDARY_CONDITION,bc);
      p.get(SOURCE,bf);

//    Instantiate equation, set terms, data and then solve equation
      LinearPDE2D eq(ms,u);
      eq.set_02();
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.setInput(SOURCE,bf);
      eq.run();

//    Get analytical solution and compute error
      Vect<double> sol(ms);
      p.get(SOLUTION,sol);
      cout << "L2-Error: " << (u-sol).Norm(WNORM2) << endl;
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}