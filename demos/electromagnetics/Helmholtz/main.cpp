/*==============================================================================

               A Finite Element Code to Solve Helmholtz Equation
                      in a Bounded Domain using OFELI

  ------------------------------------------------------------------------------

   Copyright (C) 1998 - 2024 Rachid Touzani

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

  ==============================================================================*/

#include "OFELI.h"
#include "Electromagnetics.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   void setBC(Vect<real_t>& bc);
   void error(const Vect<real_t>& u);

   ifstream mf, bcf;
   if (argc < 2) {
     cout << "\nUsage: " << argv[0] << " <parameter_file>\n";
     return EXIT_FAILURE;
   }

   IPF data("Helmholtz - 1.0",argv[1]);
   Verbosity = data.getVerbose();

   if (Verbosity) {
      cout << endl << endl;
      cout << "    *******************************************************\n";
      cout << "    *                      H e l m h o l t z              *\n";
      cout << "    *            Helmholtz equation in a bounded domain   *\n";
      cout << "    *******************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "               A Finite Element Code for Helmholtz equation\n";
      cout << "                        in a 2-D bounded domain\n\n";
      cout << "            Helmholtz uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   1.0\n\n";
      cout << "                     Copyright R. Touzani, 1999\n\n";
      cout << "=====================================================================\n\n";
   }

// ----------
// Read data
// ---------

   try {
      if (Verbosity > 3)
         cout << "Reading mesh data ...\n";

//    Read Mesh data
      Mesh ms(data.getMeshFile());
      if (Verbosity > 5)
         cout << ms;

//    Declare problem data (matrix, rhs, boundary conditions, body forces)
      if (Verbosity > 3)
         cout << "Allocating memory for matrix and R.H.S. ...\n";
      Vect<real_t> u(ms), bc(ms);

//    Read boundary conditions, body and boundary forces
      if (Verbosity > 3)
         cout << "Reading boundary conditions ...\n";
      setBC(bc);

      HelmholtzBT3 eq(ms,u);
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.set_omega(data.getString("omega"));
      eq.run();

      if (Verbosity > 4)
         cout << u;
 
      error(u);
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}


void setBC(Vect<real_t>& bc)
{
   NodeLoop(bc.getMesh()) {
      double x = TheNode.getX(), y = TheNode.getY();
      if (TheNode.getCode(1)==1) {
         bc(2*theNodeLabel-1) = cos(y)*cos(2*x);
         bc(2*theNodeLabel  ) = cos(y)*sin(2*x);
      }
   }
}


void error(const Vect<real_t>& u)
{
   real_t ee=0, e1, e2;
   size_t n = 0;
   NodeLoop(u.getMesh()) {
      n++;
      double x = TheNode.getX(), y = TheNode.getY();
      e1 = cos(y)*cos(2*x) - u(2*theNodeLabel-1);
      e2 = cos(y)*sin(2*x) - u(2*theNodeLabel);
      ee += e1*e1 + e2*e2;
   }
   cout << "Discrete L2 error = " << sqrt(ee/n) << endl;
}
