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

          Solution of the n body problem using the ODESolver class

       The n body problem can be written as a first order system of ODE:
       
           r'(i) = p(i)
           p'(i) = sum_{j!=i} G m_j (r(i)-r(j))/|r(i)-r(j)|^3

       with appropriate initial conditions

       The unknowns are stored in vector y with dimension 2*m, m = nb*dim
           y(dim*(i-1)+j) = r(i,j)
           y(m+dim*(i-1)+j) = p(i,j) 
           for 1 <= i <= nb, 1 <= j <= dim

  ==============================================================================*/

#include "OFELI.h"
#include "nBody.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc<2) {
      cout << "Usage: " << argv[0] << " <project_file>" << endl;
      return EXIT_FAILURE;
   }
   IPF data("nBody - 1.0",argv[1]);
   theTimeStep = data.getTimeStep();
   theFinalTime = data.getMaxTime();
   size_t dim = data.getInteger("dim");
   size_t nb = data.getInteger("nb");
   double G = data.getDouble("G");
   Vect<double> m(nb), y(2*dim*nb);
   m = data.getArrayDouble("m");
   y = data.getArrayDouble("init");

// Files to store trajectories and phases
   ofstream tf(data.getPlotFile(1)), pf(data.getPlotFile(2));

// Solution as a system of first-order ODEs
   try {
      nBody nbody(dim,G,m);
      ODESolver ode(nbody,RK4,theTimeStep,theFinalTime,2*dim*nb);

//    Set differential equation system
      ode.setInitial(y);

//    Loop on time steps to run the time integration scheme
      TimeLoop {

//       Run one time step
         ode.runOneTimeStep();
         tf << theTime << "    ";
         for (size_t i=1; i<=nb; ++i) {
            for (size_t j=1; j<=dim; ++j) {
               tf << y(dim*(i-1)+j) << "  ";
               pf << y(dim*(i-1)+j) << "  " << y(dim*(i-1)+j+dim*nb) << "  ";
            }
            tf << "  "; pf << " ";
         }
         tf << endl; pf << endl;
      }

   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}
