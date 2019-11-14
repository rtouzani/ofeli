/*==============================================================================

                                **********************
                                *    laplace_demo3   *
                                **********************


                  A Finite Element Code for solving the Laplace equation
                   The code uses the standard P1 finite element method

  ------------------------------------------------------------------------------

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

  ==============================================================================*/

#include "OFELI.h"
#include "Laplace.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
// Expand arguments
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <project_file>\n";
      return 0;
   }
   IPF data("laplace_demo3 - 1.0",argv[1]); 
   cout << endl << endl;
   cout << "    *******************************************************\n";
   cout << "    *               l a p l a c e - d e m o 3             *\n";
   cout << "    *   A demo program for the 3-D Laplace equation       *\n";
   cout << "    *            using P1 finite elements                 *\n";
   cout << "    *******************************************************\n\n\n";

// Read Mesh data and options
   try {
      Verbosity = data.getVerbose();
      Mesh ms(data.getMeshFile());

//    Declare problem data
      Vect<double> bc(ms), f(ms);
      Prescription p(ms,data.getDataFile());
      p.get(BOUNDARY_CONDITION,bc);
      p.get(SOURCE,f);

//    Declare problem solution and instantiate equation
      Vect<double> u(ms);
      Laplace3DT4 eq(ms,u);

//    Prescribe data
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.setInput(SOURCE,f);
      //      eq.getLinearSolver().setTolerance(1.e-8);

//    Solve problem
      eq.run();
      saveField(u,"u.pos",GMSH);

//    Get analytical solution and compute error
      Vect<double> v(ms);
      p.get(SOLUTION,v);
      cout << "L2-Error: " << (u-v).Norm(WNORM2) << endl;

   } CATCH_EXCEPTION
   return 0;
}
