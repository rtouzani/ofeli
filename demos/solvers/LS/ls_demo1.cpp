/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

  Copyright (C) 1998 - 2026 Rachid Touzani

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
 
   The tested linear system is:   [A]{x} = {b}
   where A is a tridiagonal matrix.

 ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
// Read system size
   if (argc<2) {
      cout << "Usage: " << argv[0] << " <size>" << endl;
      return EXIT_FAILURE;
   }
   size_t n=atoi(argv[1]);
   Verbosity = 4;

// Initialize matrix and right-hand side vector
   try {
      TrMatrix<double> A(n);
      Vect<double> b(n), x(n);
      for (size_t i=2; i<n; i++) {
         A(i,i)   =  2.0;
         A(i,i-1) = -1.0;
         A(i,i+1) = -1.0;
      }
      A(1,1)   =  2.0;
      A(1,2)   = -1.0;
      A(n,n-1) = -1.0;
      A(n,n)   =  2.0;
      b = 0.0;
      b(n) = n+1;

//    Solve the linear system by Gauss elimination
      LinearSolver ls(A,b,x);
      ls.setSolver(DIRECT_SOLVER);
      ls.solve();

//    Output solution vector
      cout << "Solution:\n" << x << endl; 
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}
