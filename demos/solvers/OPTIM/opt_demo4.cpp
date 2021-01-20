/*==============================================================================

                                **********************
                                *      opt_demo4     *
                                **********************


                  A demo program for solving a linear program problem 
                             using the Simplex method

  ------------------------------------------------------------------------------

   Copyright (C) 1998 - 2021 Rachid Touzani

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
 
       Maximize x1 + x2 + 3x3 -0.5x4 
       Subject to the constraints:
          x1  + 2x3 <= 740
          2x2 - 7x4 <= 0
          x2  - x3 + 2x4 >= 0.5
          x1 + x2 + x3 +x4 = 9
          and all x's >=0.

  ==============================================================================*/


#include "OFELI.h"
using namespace OFELI;

int main()
{
   int nv=4, nb_le=2, nb_ge=1, nb_eq=1;
   Vect<double> x(nv), a(nv);

   LPSolver s;
   s.setSize(nv,nb_le,nb_ge,nb_eq);
   s.set(x);
   a(1) = -1.0, a(2) = -1.0, a(3) = -3.0, a(4) = 0.5;
   s.set(LPSolver::OBJECTIVE,a);
   a(1) = 1.0, a(2) = 0.0, a(3) = 2.0, a(4) = 0.0;
   s.set(LPSolver::LE_CONSTRAINT,a,740.0);
   a(1) = 0.0, a(2) = 2.0, a(3) = 0.0, a(4) = -7.0;
   s.set(LPSolver::LE_CONSTRAINT,a,0.0);
   a(1) = 0.0, a(2) = 1.0, a(3) = -1.0, a(4) = 2.0;
   s.set(LPSolver::GE_CONSTRAINT,a,0.5);
   a(1) = 1.0, a(2) = 1.0, a(3) = 1.0, a(4) = 1.0;
   s.set(LPSolver::EQ_CONSTRAINT,a,9.0);

   s.run();
   cout << "Solution\n" << x;
   cout << "Optimal objective: " << s.getObjective() << endl;
   return 0;
}
