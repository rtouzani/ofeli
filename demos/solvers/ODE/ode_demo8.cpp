/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

  Copyright (C) 1998 - 2025 Rachid Touzani

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
       
           y'(i) = 1/m_i p(i)
           p'(i) = G (m_i m_1 )

       with appropriate initial conditions

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
   Vect<double> m(nb);
   for (size_t i=1; i<=nb; ++i)
      m(i) = data.getInteger("m"+toString(i));
   ofstream outf(data.getPlotFile(1));

// Solution as a system of first-order ODEs
   try {
      nBody nbody(dim,G,m);
      ODESolver ode(nbody,RK4,theTimeStep,theFinalTime,2*dim*nb);

//    Set differential equation system
      Vect<double> y(2*dim*nb);
      y(1) = 0.; y(2) = 0., y(3) = 1., y(4) = 0.;
      ode.setInitial(y);

//    Loop on time steps to run the time integration scheme
      TimeLoop {

//       Run one time step
         ode.runOneTimeStep();
         outf << theTime << "    ";
         for (size_t i=1; i<=nb; ++i) {
            for (size_t j=1; j<=dim; ++j)
               outf << y(dim*(i-1)+j) << "  ";
            outf << "  ";
         }
         outf << endl;
      }

   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}
