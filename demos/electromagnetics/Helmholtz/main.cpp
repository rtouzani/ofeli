/*==============================================================================

               A Finite Element Code to Solve Helmholtz Equation
                      in a Bounded Domain using OFELI

  ------------------------------------------------------------------------------

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

  ==============================================================================*/

#include "OFELI.h"
#include "Electromagnetics.h"
#include "User.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   ifstream mf, bcf;
   if (argc < 2) {
     cout << "\nUsage:  helmholtz <parameter_file>\n";
     return 0;
   }

   IPF data("Helmholtz - 1.0",argv[1]);
   int output_flag = data.getOutput();
   int bc_flag = data.getBC();
   double wave_nb = data.getDoublePar(1);

   if (output_flag) {
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

// Read Mesh data
   try {
      if (output_flag > 1)
         cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile());
      wave_nb = data.getDouble("wave_nb");

      if (output_flag > 1)
         cout << ms;
      User ud(ms);

//    Declare problem data (matrix, rhs, boundary conditions, body forces)
      if (output_flag > 1)
         cout << "Allocating memory for matrix and R.H.S. ...\n";
      SkMatrix<complex_t> A(ms);
      Vect<complex_t> b(ms.getNbDOF());

//    Read boundary conditions, body and boundary forces
      if (output_flag > 1)
         cout << "Reading boundary conditions ...\n";
      Vect<complex_t> bc(ms.getNbDOF());
      if (!bc_flag)
         ud.setDBC(bc);

//    Read in boundary conditions and body forces
      ud.setDBC(bc);

//    Loop over elements
//    ------------------

      MeshElements(ms) {
         HelmholtzBT3 eq(theElement);
         eq.LHS(wave_nb);
         eq.ElementAssembly(A);
      }

//    Loop over sides
//    ---------------

      MeshSides(ms) {
         HelmholtzBT3 eq(theSide);
         eq.BoundaryRHS(ud);
         eq.SideAssembly(b);
      }

      A.Prescribe(b,bc);
      A.solve(b);

      if (output_flag > 0)
         cout << b;

      void error(const Mesh &ms, User &ud, const Vect<complex<double> > &u);
      error(ms,ud,b);
   } CATCH_EXCEPTION
   return 0;
}
