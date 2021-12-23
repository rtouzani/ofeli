/*==============================================================================

                              **********************
                              *    T  I  F  F  2   *
                              **********************

                 A Finite Element Code for Transient Incompressible
                       Fluid Flow Simulations in 2-D Geometries

  ------------------------------------------------------------------------------

   Copyright (C) 1998 - 2022 Rachid Touzani

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
#include "Fluid.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <parameter_file>\n";
      return 0;
   }

   IPF data("tiff2 - 1.0",argv[1]);
   Verbosity = data.getVerbose();
   int save_flag = data.getSave();
   theTimeStep = data.getTimeStep();
   theFinalTime = data.getMaxTime();

   if (Verbosity) {
     cout << "=====================================================================\n\n";
     cout << "                             T  I  F  F  2\n\n";
     cout << "         A Finite Element Code for Transient Incompressible\n";
     cout << "                 Fluid Flow Simulation in 2-D Geometries\n\n\n";
     cout << "              tiff2 uses OFELI Library of Finite Element Classes\n\n";
     cout << "                             V E R S I O N   1.0\n\n";
     cout << "                       Copyright R. Touzani, 1998\n\n";
     cout << "=====================================================================\n\n";
   }

//-----------
// Read data
//-----------

// Read Mesh data
   try {
      if (Verbosity > 1)
        cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile());

//    Declare problem data (boundary conditions, body forces)
      Vect<double> u(ms), bc(ms), p(ms,NODE_DOF,1);
      cout << "**** "<<p.getNbDOF()<<endl;
      Prescription pr(ms,data.getDataFile());
      pr.get(BOUNDARY_CONDITION,bc);
      u.setName("Velocity");

      IOField v_file(data.getMeshFile(),data.getString("v_file"),ms,IOField::OUT);
      IOField p_file(data.getMeshFile(),data.getString("p_file"),ms,IOField::OUT);

//    Loop over time steps
//    --------------------

      NSP2DQ41 eq(ms,u);
      eq.setInput(BOUNDARY_CONDITION,bc);
      eq.setInput(PRESSURE_FIELD,p);

      TimeLoop {

         if (Verbosity > 1)
            cout << "Performing time step " << theStep << " ..." << endl;
         cout << "Step: " << theStep << ", Time: " << theTime << endl;

//       Solve for the present time step
         eq.runOneTimeStep();

//       Store solution
         if (save_flag) {
            v_file.put(u);
            p_file.put(p);
	 }
      }
   } CATCH_EXCEPTION
   return 0;
}
