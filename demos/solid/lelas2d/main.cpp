/*==============================================================================

                             *************************
                             *     L E L A S 2 D     *
                             *************************

                A Finite Element Code for Linearized Elastostatics in
                                  Two Dimensions

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
#include "Solid.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <parameter_file>\n";
      return 0;
   }

   IPF data("lelas2d - 1.1 ",argv[1]);
   int output_flag = data.getOutput();
   int save_flag = data.getSave();

   if (output_flag) {
      cout << endl << endl;
      cout << "    *******************************************************\n";
      cout << "    *                   L E L A S 2 D                     *\n";
      cout << "    *     Two-Dimensional Linearized Elastostatics        *\n";
      cout << "    *******************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "         A Finite Element Code for Linearized Elastostatics\n";
      cout << "                             in 2-D Geometries\n\n";
      cout << "            LELAS2D uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   1.0\n\n";
      cout << "                       Copyright R. Touzani, 1998\n\n";
      cout << "=====================================================================\n\n";
   }

//----------
// Read data
//----------

// Read Mesh data
   try {
      if (output_flag > 1)
         cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile(1));
      Prescription p(ms,data.getDataFile());
      if (output_flag > 1)
         cout << ms;

//    Read boundary conditions, body and boundary forces
      Vect<double> u(ms), bc(ms), body_f(ms), bound_f(ms,2,SIDE_DOF);
      p.get(BOUNDARY_CONDITION,bc,0);
      if (output_flag > 1)
         cout << "Reading body forces ..." << endl;
      p.get(BODY_FORCE,body_f);
      if (output_flag > 1)
         cout << "Reading Boundary Tractions ..." << endl;
      p.get(TRACTION,bound_f,0);

//    Declare equation instance and solve
      if (output_flag > 1)
         cout << "Setting and solving equation ...\n";
      Elas2DQ4 eq(ms,u);
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.setInput(BODY_FORCE,body_f);
      eq.setInput(TRACTION,bound_f);
      eq.run();

//    Output and save solution
      if (output_flag > 0)
         cout << u;
      if (save_flag) {
         IOField pl_file(data.getPlotFile(),IOField::OUT);
         pl_file.put(u);
      }

//    Calculate principal and Von-Mises stresses
      if (output_flag > 1)
         cout << "Calculating stresses ...\n";
      Vect<double> st, vm;
      eq.Stress(st,vm);
   } CATCH_EXCEPTION
   return 0;
}
