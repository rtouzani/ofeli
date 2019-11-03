/*==============================================================================

                                *******************
                                *     S T D 3     *
                                *******************


                        A Finite Element Code for Steady-State
                        Analysis of Thermal Diffusion Problems
                                  in 3-D Geometries


  ------------------------------------------------------------------------------

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

  ==============================================================================*/


#include "OFELI.h"
#include "Therm.h"
using namespace OFELI;


int main(int argc, char *argv[])
{
// Expand arguments
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <param_file>\n";
      return 0;
   }

   IPF data("std3 - 1.0",argv[1]);
   Verbosity = data.getVerbose();
   int save_flag = data.getSave();

   if (Verbosity) {
      cout << endl << endl;
      cout << "    *******************************************************\n";
      cout << "    *                      S  T  D  3                     *\n";
      cout << "    *             Steady State Thermal Diffusion          *\n";
      cout << "    *                   in 3-D Geometries                 *\n";
      cout << "    *******************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "               A Finite Element Code for Steady-State\n";
      cout << "                Analysis of Thermal Diffusion Problems\n";
      cout << "                             in 3-D Geometries\n\n";
      cout << "            STD3 uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   1.0\n\n";
      cout << "                     Copyright R. Touzani, 2002\n\n";
      cout << "=====================================================================\n\n";
   }

//---------------------------------
// Read data
//---------------------------------

   try {

// Read Mesh data
      if (Verbosity>2)
         cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile());
      Prescription p(ms,data.getDataFile());
      if (Verbosity>3)
         cout << ms;

//    Declare problem data (matrix, rhs, boundary conditions, body forces)
      Vect<double> u(ms);

//    Read boundary conditions, body and boundary forces
      if (Verbosity>1)
         cout << "Reading boundary conditions, body and boundary sources ...\n";
      Vect<double> bc(ms), bf(ms), sf(ms,1,BOUNDARY_SIDE_FIELD);
      p.get(BOUNDARY_CONDITION,bc);
      p.get(BODY_FORCE,bf);
      p.get(BOUNDARY_FORCE,sf,0);

//    Instantiating equation
      if (Verbosity>1)
         cout << "Setting equation" << endl;
      DC3DT4 eq(ms,u);
      if (Verbosity>1)
         cout << "Setting equation data ..." << endl;
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.setInput(SOURCE,bf);
      eq.setInput(FLUX,sf);
      eq.setTerms(DIFFUSION);
      eq.setSolver(CG_SOLVER,DILU_PREC);

//    Run solution
      if (Verbosity>1)
         cout << "Solving equation" << endl;
      eq.run();
      if (Verbosity>4)
         cout << u;

      if (save_flag) {
         if (Verbosity>1)
            cout << "Saving solution" << endl;
         saveField(u,"beam.pos",GMSH);
      }
   } CATCH_EXCEPTION
   return 0;
}
