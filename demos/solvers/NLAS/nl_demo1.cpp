/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

  Copyright (C) 1998 - 2023 Rachid Touzani

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
    The problem is solved by the Regula Falsi method.

    In the example f(x) = x^2-2.
    The function f is given as a user defined function

 ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;

double f(double x)
{
   return x*x-2;
}

int main(int argc, char *argv[])
{
   double x = 0.;
   NLASSolver nls(x,REGULA_FALSI);
   nls.setFunction(function<double(double)>(f));
   Verbosity = 3;
   nls.setInitial(0.,2.);
   nls.run();
   cout << "Solution: " << x << endl;
   cout << nls;
   return EXIT_SUCCESS;
}
