/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

  Copyright (C) 1998 - 2018 Rachid Touzani

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
                     Case of a multi-variable problem

      The objective function is 
      f(x1,x2,x3) = x1^2 + (x1-x2)^2 + (x2-x3)^2 + x3^2 - 8*x3

      The minimum is achieved for x1 = 1, x2 = 2, x3 = 3

 ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;


int main()
{
   try {
//    Instantiate solution vector (Initialized to 0)
      Vect<double> x(3);

//    Instantiate optimization solver class using solution vector
      OptSolver os(x);

//    Select optimization algorithm
      os.setOptMethod(OptSolver::GRADIENT);

//    Set objective function and its gradient 
      os.setObjective("x1^2+(x1-x2)^2+(x2-x3)^2+x3^2-8*x3");
      os.setGradient("4*x1-2*x2",1);
      os.setGradient("-2*x1+4*x2-2*x3",2);
      os.setGradient("-2*x2+4*x3-8",3);

//    Run the optimization procedure
      os.run();

//    Output class information and solution
      cout << os;
      cout << "\nSolution:\n" << x;
   } CATCH_EXCEPTION

   return 0;
}
