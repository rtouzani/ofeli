/*==============================================================================

                                **********************
                                *    laplace_demo2   *
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
   IPF data("laplace_demo2 - 1.0",argv[1]); 
   cout << endl << endl;
   cout << "    *******************************************************\n";
   cout << "    *               l a p l a c e - d e m o 2             *\n";
   cout << "    *       A demo program for the Laplace equation       *\n";
   cout << "    *            using P2 finite elements                 *\n";
   cout << "    *******************************************************\n\n\n";

// Read Mesh data and options
   try {
      Verbosity = data.getVerbose();
      Mesh ms(data.getMeshFile());

// Transform 3-Node to 6-Node triangles and save resulting mesh
      ms.AddMidNodes();
      ms.put("mesh-P2.m");

//    Declare problem data
      Vect<double> bc(ms), bf(ms);
      Prescription p(ms,data.getPrescriptionFile());
      p.get(BOUNDARY_CONDITION,bc);
      p.get(SOURCE,bf);


//    Declare problem solution and instantiate equation
      Vect<double> u(ms);
      Laplace2DT6 eq(ms,u);

//    Prescribe data
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.setInput(SOURCE,bf);

//    Solve problem
      eq.run();

//    Get analytical solution and compute error
      Vect<double> sol(ms);
      p.get(SOLUTION,sol);
      cout << "L2-Error: " << (u-sol).Norm(WNORM2) << endl;

   } CATCH_EXCEPTION
   return 0;
}
