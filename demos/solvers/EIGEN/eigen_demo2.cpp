/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2026 Rachid Touzani

   This program is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free
   Software Foundation; Version 2 of the License.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the :

   Free Software Foundation
   Inc., 59 Temple Place - Suite 330
   Boston, MA  02111-1307, USA

  ==============================================================================

                       An example of the eigenproblem solver

      Solution of an eigenvalue problem for a given dense unsymmetric matrix

  ==============================================================================*/

#include "OFELI.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   banner("eigen_demo2");
   const int n = 4;

// Define the matrix
   try {
      DMatrix<double> A(n,n);
      A = 0;
      A(1,1) =  2.0; A(1,2) = -1.0;
      A(2,1) = -1.0; A(2,2) =  2.0; A(2,3) = -1.0;
      A(3,2) = -1.0; A(3,3) =  2.0; A(3,4) = -1.0;
      A(4,1) =  1.0; A(4,2) =  1.0; A(4,3) =  1.0; A(4,4) = 4.0;

//    Solve the eigenvalue problem
      Vect<double> evr, evi;
      EigenProblemSolver e(A,true);

//    Output eigenvalues and eigenvectors
      Verbosity = 10;
      cout << e;
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}
