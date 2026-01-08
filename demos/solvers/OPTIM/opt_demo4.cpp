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

              Solution of the well known Brachistochrone problem
  
  ==============================================================================*/

#include "OFELI.h"
#include "Opt4.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
// Define number of discretization points and bound values
   int n = 20;
   const double a = 0., b = 1., g = 10., y0=5., y1=1.;
   if (argc>1)
      n = atoi(argv[1]);

   try {
//    Instantiate solution vector (Initialized to 0)
      Vect<double> y(n+1);
      y[0] = y0, y[n] = y1;

//    Instantiate optimization solver class using solution vector
      Verbosity = 1;
      Opt4 opt(g,a,b,y);
      OptSolver os(opt,y);

//    Select optimization algorithm
      os.setOptMethod(OptSolver::TRUNCATED_NEWTON);

//    Impose bound constraints
      os.setEqBound(1,y0);
      os.setEqBound(n+1,y1);

//    Run the optimization procedure
      os.run();

//    Output class information, solution and save for plotting
      cout << os;
      cout << "\nSolution stored in file 'plot.dat'" << endl;
      double h = (b-a)/n;
      ofstream fs("plot.dat");
      for (int i=0; i<=n; ++i) 
         fs << i*h << "  " << y[i] << endl;
   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}
