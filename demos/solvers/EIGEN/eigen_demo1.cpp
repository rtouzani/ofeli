/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

       Solution of an eigenvalue problem for a given dense symmetric matrix

  ==============================================================================*/

#include "OFELI.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   banner("eigen_demo1");
   const int n = 6;

// Define the matrix
   try {
      DSMatrix<double> A(n);
      A = 0;
      A(1,1) = 2; A(1,2) = -1;
      for (size_t i=2; i<n; i++) {
         A(i,i-1) = -1;
         A(i,i) = 2;
         A(i,i+1) = -1;
      }
      A(n,n-1) = -1; A(n,n) = 2;

//    Solve the eigenvalue problem
      Vect<double> ev;
      EigenProblemSolver e(A,ev);

//    Output eigenvalues, save eigenvectors
      Vect<double> v(n);
      for (int i=1; i<=n; i++) {
         cout << "Eigenvalue #" << i << ": " << ev(i) << endl;
         e.getEigenVector(i,v);
         cout << "Eigen vector:\n" << v;
      }
      cout << "Nb. of iterations: " << e.getNbIter() << endl;
   } CATCH_EXCEPTION
   return 0;
}
