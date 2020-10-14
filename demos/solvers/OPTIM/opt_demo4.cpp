/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                A program to illustrate the Optimization solver

        Solution of a Linear Programming problem using the Simplex method
 
   A simple production planning problem is given by the use of two ingredients 
   A and B that produce products 1 and 2. The available supply of A is 30 units
   and B is 44 units. For production it requires:

       3 units of A and 8 units of B to produce Product 1
       6 units of A and 4 units of B to produce Product 2 

   There are at most 5 units of Product 1 and 4 units of Product 2. Product 1 
   can be sold for 100 and Product 2 can be sold for 125. The objective is to 
   maximize the profit for this production problem. 

   The linear program can be stated as follows:

   Maximize profit : 100*x1 + 125*x2
   Subject to the constraints:
        3*x1 + 6*x2 <= 30, 8*x1 + 4*x2 <= 44, 0<x1<5, 0<x2<4

  ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   try {
#ifdef EX1
    const size_t nc = 6, nr = 3;
    DMatrix<double> A(nr,nc);
    Vect<double> b(nr), c(nc), x(nr);
    c(1) = -6., c(2) = -5, c(3) = -4, c(4) = 0., c(5) = 0., c(6) = 0.;
    A(1,1) = 2., A(1,2) = 1., A(1,3) = 1., A(1,4) = 1., A(1,5) = 0., A(1,6) = 0.; 
    A(2,1) = 1., A(2,2) = 3., A(2,3) = 2., A(2,4) = 0., A(2,5) = 1., A(2,6) = 0.; 
    A(3,1) = 2., A(3,2) = 1., A(3,3) = 2., A(3,4) = 0., A(3,5) = 0., A(3,6) = 1.; 
    b(1) = 180., b(2) = 300., b(3) = 240.;
#endif
#ifdef EX2
    const size_t nc = 2, nr = 2;
    DMatrix<double> A(nr,nc);
    Vect<double> b(nr), c(nc), x(nr);
    c(1) = -5., c(2) = -7.;
    A(1,1) = 3., A(1,2) = 4.; 
    A(2,1) = 2., A(2,2) = 3.;
    b(1) = 650., b(2) = 500.;
#endif
#ifdef EX3
    const size_t nc = 6, nr = 4;
    DMatrix<double> A(nr,nc);
    Vect<double> b(nr), c(nc), x(nc);
    c(1) = -2., c(2) = -1., c(3) = 0., c(4) = 0., c(5) = 0., c(6) = 0.;
    A(1,1) = -1., A(1,2) = -1., A(1,3) = 1., A(1,4) = 0., A(1,5) = 0., A(1,6) = 0.; 
    A(2,1) = -1., A(2,2) =  1., A(2,3) = 0., A(2,4) = 1., A(2,5) = 0., A(2,6) = 0.;
    A(3,1) = -1., A(3,2) =  0., A(3,3) = 0., A(3,4) = 0., A(3,5) = 1., A(3,6) = 0.;
    A(4,1) =  0., A(4,2) = -1., A(4,3) = 0., A(4,4) = 0., A(4,5) = 0., A(4,6) = 1.;
    b(1) = -6., b(2) = -4., b(3) = 0., b(4) = 0.;
#endif

      Simplex s(A,b,c,x);
      s.run();
      cout << "Solution:\n" << x;
      cout << "Cost: " << s.getMaxCost() << endl;
      cout << "Nb. of iterations: " << s.getNbIter() << endl;
   } CATCH_EXCEPTION

   return 0;
}
