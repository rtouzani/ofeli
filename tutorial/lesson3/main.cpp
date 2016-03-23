/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

           Solution of a 2-D Diffusion problem using P1 Finite elements
           Linear System is solved using the Conjugate Gradient Method

  This function illustrates the use of an iterative solver for three
  distinct situations:

      1. Using native iterative solvers of OFELI
      2. Using eigen library
      3. Using Petsc library
  ==============================================================================*/

#include "OFELI.h"
#include "Therm.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   banner();

// Read and output mesh data
   if (argc<=1) {
      cout << "Usage: lesson3 <mesh_file>" << endl;
      exit(1);
   }
   Mesh ms(argv[1],true);

// Declare problem data (matrix, rhs, boundary conditions, body forces)
#if defined(USE_PETSC)
   PETScWrapper<double> w(1,argv);
   PETScMatrix<double> A(ms);
   PETScVect<double> b(ms.getNbEq()), x(ms.getNbEq()), bc(ms), u(ms);
#else
   SpMatrix<double> A(ms);
   Vect<double> b(ms.getNbEq()), x(ms.getNbEq()), bc(ms), u(ms);
#endif
   bc.setNodeBC(1,"y");

// Loop over elements
// ------------------

   MeshElements(ms) {

//    Declare an instance of class DC2DT3
      DC2DT3 eq(theElement);

//    Diffusion contribution to matrix
      eq.Diffusion();

//    Boundary condition contribution to RHS
      eq.updateBC(bc);

//    Assemble element matrix and RHS
      eq.ElementAssembly(A);
      eq.ElementAssembly(b);
   }

// Solve the linear system of equations by an iterative method
#if defined(USE_EIGEN)
   LinearSolver<double> ls(1000,1.e-8,1);
   int nb_it = ls.solve(A,b,x,CG_SOLVER,DIAG_PREC);
#elif defined(USE_PETSC)
   w.setLinearSystem(A,b,KSPCG,PCJACOBI,1.e-8);
   w.solve(x);
   int nb_it = w.getIterationNumber();
#else
   LinearSolver<double> ls(1000,1.e-8,0);
   int nb_it = ls.solve(A,b,x,CG_SOLVER,DILU_PREC);
#endif
   cout << "Number of iterations: " << nb_it << endl;

// Output solution
   u.insertBC(ms,x,bc);
   return 0;
}
