/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2016 Rachid Touzani

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

               An example of a Finite Element Code using OFELI

            Solution of a 2-D Poisson problem using P1 Finite elements
  A truncated Newton Algorithm is used to solve the resulting optimization problem

  ==============================================================================*/


#include "OFELI.h"
#include "Opt.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc <= 1) {
      cout << " Usage : lesson5 <par_file>\n";
      exit(1);
   }
   IPF data(argv[1]);
   Mesh ms(data.getMeshFile(1));
   int n = ms.getNbDOF();

// Declare solution vector
   Vect<double> x(n),low(n),up(n);
   Vect<int> pivot(n);

// Read in boundary conditions
   Vect<double> bc(n);
   bc.setNodeBC(ms,2,"1");

// Define and Solve the Optimization Problem
   Opt theOpt(ms);
   x = 0;
   BCAsConstraint(ms,bc,up,low);
   OptimTN<Opt>(theOpt,x,low,up,pivot,100,1.e-12,1);

// Output solution
   cout << "\nSolution:\n" << x;

   return 0;
}
