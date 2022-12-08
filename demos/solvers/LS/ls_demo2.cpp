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

                A program to illustrate the LinearSolver class
                           for solving linear systems
 
   The tested linear system issues from a finite element discretization of the
   Poisson problem.
   Since the resulting matrix is symmetric positive definite, we use the
   Conjugate Gradient method preconditioned by the Diagonal preconditioner

 ==============================================================================*/

#include "OFELI.h"
#include "Laplace.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
// Read and output mesh data
   if (argc<=1) {
      cout << "Usage: " << argv[0] << " <mesh_file>" << endl;
      return EXIT_FAILURE;
   }
   Verbosity = 4;

   try {
      Mesh ms(argv[1],true);
      Vect<double> u(ms), bc(ms), f(ms);
      f = "exp(-20*(x*x+y*y))";
      bc = 0;

      Laplace2DT3 eq(ms);
      eq.setInput(SOLUTION,u);
      eq.setInput(SOURCE,f);
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.setSolver(CG_SOLVER,DILU_PREC);
      eq.run();
      cout<<u;
      saveField(u,"sol.pos",GMSH);
      cout << "Solution is stored in file 'sol.pos'" << endl;
      cout << "You can plot this using GMSH" << endl;
   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}
