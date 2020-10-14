/*==============================================================================

                                **********************
                                *    lp_demo1   *
                                **********************


                  A demo program for solving a linear program problem 
                             using the Simplex method

  ------------------------------------------------------------------------------

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

  A simple production planning problem is given by the use of two ingredients A 
  and B that produce products 1 and 2. The available supply of A is 30 units 
  and B is 44 units. For production it requires:

    3 units of A and 8 units of B to produce Product 1
    6 units of A and 4 units of B to produce Product 2 

  There are at most 5 units of Product 1 and 4 units of Product 2. Product 1 can 
  be sold for 100 and Product 2 can be sold for 125. The objective is to maximize 
  the profit for this production problem. 

  The linear program can be stated as follows:
      Maximize profit : 100*x1 + 125*x2
      Subject to the constraints:
      3*x1 + 6*x2 <= 30, 8*x1 + 4*x2 <= 44, 0<x1<5, 0<x2<4

  ==============================================================================*/


#include "OFELI.h"
using namespace OFELI;


int main()
{
   const size_t n = 2, m = 6;
   DMatrix<double> A(m,n);
   Vect<double> c(n), b(m), x(n+m);
   c(1) = -100., c(2) = -125.;
   b(1) = 30., b(2) = 44., b(3) = 0., b(4) = 5., b(5) = 0., b(6) = 4.;
   A(1,1) = 3., A(1,2) = 6., A(2,1) = 8., A(2,2) =  4., A(3,1) = -1., A(3,2) = 0.;
   A(4,1) = 1., A(4,2) = 0., A(5,1) = 0., A(5,2) = -1., A(6,1) =  0., A(6,2) = 1.;
   Simplex s(A,b,c,x);
   int nb_it = s.run();
   cout << "Solution:\n" << x;
   cout << "Cost: " << s.getCost() << endl;
   cout << "Nb. of iterations: " << nb_it << endl;
   return 0;
}
