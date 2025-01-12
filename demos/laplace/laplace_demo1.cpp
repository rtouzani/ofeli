/*==============================================================================

                                **********************
                                *    laplace_demo1   *
                                **********************


                  A Finite Element Code for solving the Laplace equation
                   The code uses the standard P1 finite element method

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
   IPF data("laplace_demo1 - 1.0",argv[1]); 
   cout << endl << endl;
   cout << "    *******************************************************\n";
   cout << "    *               l a p l a c e - d e m o 1             *\n";
   cout << "    *       A demo program for the Laplace equation       *\n";
   cout << "    *              using P1 finite elements               *\n";
   cout << "    *******************************************************\n\n\n";

   try {

//    Read Mesh data and options
      Verbosity = data.getVerbose();
      Mesh ms(data.getMeshFile());

//    Declare problem data
      Vect<double> bc(ms), bf(ms), sf(ms);
      Prescription p(ms,data.getPrescriptionFile());
      p.getBoundaryCondition(bc);
      p.getBodyForce(bf);
      p.getFlux(sf);

//    Declare problem solution and instantiate equation
      Vect<double> u(ms);
      Laplace2DT3 eq(ms,u);

//    Prescribe data
      eq.setBoundaryCondition(bc);
      eq.setSource(bf);
      eq.setSolver(CG_SOLVER,SSOR_PREC);
      eq.LinearSystemInfo();

//    Solve problem
      eq.run();
      saveField(u,"u.pos",GMSH);

//    Get analytical solution and compute error
      Vect<double> sol(ms);
      p.getSolution(sol);
      cout << "L2-Error: " << (u-sol).Norm(WNORM2) << endl;

   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}
