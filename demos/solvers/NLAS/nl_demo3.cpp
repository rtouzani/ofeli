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

                    A program to illustrate the NLAS solver

                 Solution of a system of 2 nonlinear equations
                       by the Newton's method

         x1*x1-2*x1*x2 = 2
         x1 + x2^2 = -1

     The approximate solutions the system are x1 = -1.11509, x2 =  0.339246 
                                          and x1 = -3.93432, x2 = -1.71298

 ==============================================================================*/

#include "OFELI.h"
using namespace OFELI;


int main(int argc, char *argv[])
{
   try {
      Vect<double> x(2);

//    First solution using a first initial guess
      x(1) = -1.; x(2) = 0.;
      NLASSolver nls(x,NEWTON);
      nls.setf("x1*x1-2*x1*x2-2");
      nls.setDf("2*x1-2*x2",1,1);
      nls.setDf("-2*x1",1,2);
      nls.setf("x1+x2*x2+1");
      nls.setDf("1",2,1);
      nls.setDf("2*x2",2,2);
      nls.run();
      cout << "Solution 1:\n" << x << endl;

//    Second solution using a second initial guess
      x(1) = -3.; x(2) = -1.;
      nls.setInitial(x);
      nls.run();
      cout << "Solution 2:\n" << x << endl;

   } CATCH_EXCEPTION
   return 0;
}
