/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

  Copyright (C) 1998 - 2017 Rachid Touzani

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
      cout << "Usage: ls_demo2 <mesh_file>" << endl;
      exit(1);
   }
   Mesh ms(argv[1],true);

   SpMatrix<double> A(ms);
   Vect<double> b(ms.getNbEq()), x(ms.getNbEq()), bc(ms), f(ms);
   f = "exp(-20*(x*x+y*y))";
   bc = 0;

// Loop over elements
   MeshElements(ms) {
      Laplace2DT3 eq(theElement);
      eq.LHS();
      eq.BodyRHS(f);
      eq.updateBC(bc);
      eq.ElementAssembly(A);
      eq.ElementAssembly(b);
   }

// Solve the linear system of equations by an iterative method
   LinearSolver<double> ls(1000,1.e-8,0);
   int nb_it = ls.solve(A,b,x,CG_SOLVER,DIAG_PREC);
   cout << "Number of iterations: " << nb_it << endl;

// Create solution vector by inserting boundary nodes (Vector u)
   Vect<double> u(ms);
   u.insertBC(ms,x,bc);
   return 0;
}
