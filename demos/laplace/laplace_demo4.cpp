/*==============================================================================

                                **********************
                                *    laplace_demo4   *
                                **********************


           A Boundary Element Code for solving the Steklov-Poincare problem
                The code uses the standard P0 boundary element method

  ------------------------------------------------------------------------------

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

  ==============================================================================*/

#include "OFELI.h"
#include "Laplace.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
// Expand arguments
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <project_file>\n";
      return EXIT_FAILURE;
   }
   IPF data("laplace_demo4 - 1.0",argv[1]); 
   cout << endl << endl;
   cout << "    *******************************************************\n";
   cout << "    *               l a p l a c e - d e m o 4             *\n";
   cout << "    *   A demo program for the Steklov-Poincare problem   *\n";
   cout << "    *              using P0 boundary elements             *\n";
   cout << "    *******************************************************\n\n\n";

// Read Mesh data and options
   try {
      Verbosity = data.getVerbose();
      Mesh ms(data.getMeshFile());

//    Declare problem data
      Vect<double> g(ms), u(ms,BOUNDARY_SIDE_DOF);
      g.setNodeBC(1,"0.");
      g.setNodeBC(2,"1.");

//    Declare problem solution and instantiate equation
      SteklovPoincare2DBE eq(ms,u);

//    Prescribe normal derivative on boundary as data
      eq.setBoundaryCondition(g);

//    Solve problem
      eq.run();

//    Output solution norm
      cout << "L2 solution norm: " << u.Norm(WNORM2) << endl;

   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}