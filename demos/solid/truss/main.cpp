/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2020 Rachid Touzani

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

                  Solution of a structural truss problem

  ==============================================================================*/

#include "OFELI.h"
#include "Solid.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <project_file>\n";
      return 0;
   }
   IPF data("truss",argv[1]);
   int output_flag = data.getOutput();

   if (output_flag) {
      cout << endl << endl;
      cout << "    *******************************************************\n";
      cout << "    *                      T R U S S                      *\n";
      cout << "    *   Two-Dimensional Trusses in Structural Mechanics   *\n";
      cout << "    *******************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "         A Finite Element Code for planar trusses of bars\n";
      cout << "              beam uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   1.0\n\n";
      cout << "                       Copyright R. Touzani, 2008\n\n";
      cout << "=====================================================================\n\n";
   }

// Read Mesh data
   try {
      if (output_flag > 1)
         cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile());
      if (output_flag > 1)
         cout << ms;
      Prescription p(ms,data.getDataFile());

//    Declare problem data (boundary conditions, loads)
      if (Verbosity > 1)
         cout << "Reading boundary conditions and loads ...\n";
      Vect<double> u(ms), bc(ms), f(ms);
      p.get(BOUNDARY_CONDITION,bc);
      p.get(POINT_FORCE,f);
      double section = data.getDouble("section");

//    Solve problem
      if (Verbosity > 1)
         cout << "Calculating and assembling element arrays ..." << endl;
      Bar2DL2 eq(ms,u);
      eq.setSection(section);
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.setInput(POINT_FORCE,f);
      eq.run(); 

//    Output solution
      cout << "\nSolution:\n" << u;

//    Transform mesh to deformed one
      IOField pf(data.getMeshFile(),data.getPlotFile(),ms,IOField::OUT);
      pf.put(u);
      DeformMesh(ms,u,40);
      ms.put(data.getProject()+"-1.m");
   } CATCH_EXCEPTION
   return 0;
}
