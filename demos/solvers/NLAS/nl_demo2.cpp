/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

  Copyright (C) 1998 - 2019 Rachid Touzani

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

                    A program to illustrate the NLAS solver

                      Solution of the equation    f(x) = 0
    where x is a real variable.
    The problem is solved by the secant method (Newton's method with constant
    derivative).

    In the example f(x) = x^2-2.
    The function f is given by a regular expression

 ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;


int main(int argc, char *argv[])
{
   double x = 1.;
   NLASSolver nls(x,SECANT);
   nls.setf("x*x-2");
   nls.setDf("2*x");
   nls.run();
   cout << "Solution: " << x << endl;
   cout << nls;
   return 0;
}
