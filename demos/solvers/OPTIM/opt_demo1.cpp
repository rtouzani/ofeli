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

               A program to illustrate the Optimization solver

         Case of a one variable problem with the objective function given
         by an algebraic expression

      The objective function is f(x) = x^2 - x + 1
      The minimum is achieved for x = 0.5

 ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;


int main(int argc, char *argv[])
{
// Instantiate optimization solver class
   try {
      OptSolver os;

//    Select optimization algorithm
//    Here we choose the truncated Newton method
      os.setOptMethod(OptSolver::GRADIENT);

//    Define objective function
      os.setObjective("x^2-x+1");

//    Run the optimization procedure
      os.run();

//    Output class information and solution
      cout << os << endl;
      cout << "Solution: " << os.getSolution() << endl;
   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}
