/*==============================================================================

                             *************************
                             *     C O N T A C T     *
                             *************************


                            A Signorini Contact Problem
                      for Plane Strain Linearized Elasticity

  ------------------------------------------------------------------------------

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

  ==============================================================================*/

#include "OFELI.h"
#include "Solid.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   ifstream mf, bcf, bodyf, boundf;
   ofstream dmf;

   cout << "\n\n";
   cout << "contact, version 1.1, Copyright (c) 1998 - 2008 by Rachid Touzani\n";
   cout << "contact comes with ABSOLUTELY NO WARRANTY.\n";
   cout << "This is free software, and your are allowed to redistribute it\n";
   cout << "under certain conditions. Details are distributed with the software." << endl;

   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <parameter_file>\n";
      return EXIT_FAILURE;
   }
   IPF data("contact - 1.0",argv[1]);
   Verbosity = data.getVerbose();

   if (Verbosity) {
      cout << endl << endl;
      cout << "    *********************************************************\n";
      cout << "    *                     C o n t a c t                     *\n";
      cout << "    * Two-Dimensional Linearized Elastostatics with contact *\n";
      cout << "    *********************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "         A Finite Element Code for Linearized Elastostatics\n";
      cout << "                      in 2-D Geometries with contact\n\n";
      cout << "          Contact uses OFELI Library of Finite Element Classes\n\n";
      cout << "                         V E R S I O N   1.1\n\n";
      cout << "                     Copyright R. Touzani, 2003\n\n";
      cout << "=====================================================================\n\n";
   }

//----------
// Read data
//----------

// Read parameters and mesh data
   try {
      Mesh ms(data.getMeshFile(1));
      Prescription p(ms,data.getDataFile());
      if (Verbosity > 1)
         cout << "Reading mesh data ...\n";
      if (Verbosity > 4)
         cout << ms;

//    Declare problem data (matrix, rhs, boundary conditions, body forces)

//    Read boundary conditions, body and boundary forces
      if (Verbosity > 1)
         cout << "Reading boundary conditions, body and boundary forces ...\n";
      Vect<double> u(ms), sf(ms,SIDE_DOF);
      p.getTraction(sf);

      Elas2DT3 eq(ms,u);
      eq.setBoundaryCondition(p.getBoundaryCondition());
      eq.setTraction(sf);
      eq.setBodyForce(p.getBodyForce());
      eq.run();

//    Compute deformed mesh
      ms.Deform(u);
      ms.put(data.getProject()+"-1.m");

//    Output and/or save solution
      if (Verbosity > 4)
         cout << u;
      if (data.getSave())
         saveField(u,data.getProject()+".pos",GMSH);
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}
