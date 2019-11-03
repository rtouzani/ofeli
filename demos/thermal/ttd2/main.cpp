/*==============================================================================

                                *******************
                                *     T T D 2     *
                                *******************


                        A Finite Element Code for Transient
                       Analysis of Thermal Diffusion Problems
                                  in 2-D Geometries

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
   ifstream mf, inf, bcf, bodyf, boundf;

// Expand arguments
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <parameter_file>\n";
      return 0;
   }

   IPF data("ttd2 - 1.2",argv[1]);
   Verbosity = data.getVerbose();
   theFinalTime = data.getMaxTime();
   theTimeStep = data.getTimeStep();
   int save_flag = data.getSave();
   IOField pf(data.getPlotFile(),IOField::OUT);

   cout << endl << endl;
   cout << "    *******************************************************\n";
   cout << "    *                        T  T  D  2                   *\n";
   cout << "    *    Transient Thermal Diffusion in 2-D Geometries    *\n";
   cout << "    *******************************************************\n\n\n";
   cout << "=====================================================================\n\n";
   cout << "                   A Finite Element Code for Transient\n";
   cout << "       Analysis of Thermal Diffusion Problems in 2-D Geometries\n\n";
   cout << "            ttd2 uses OFELI Library of Finite Element Classes\n\n";
   cout << "                           V E R S I O N   1.0\n\n";
   cout << "                     Copyright R. Touzani, 1998\n\n";
   cout << "=====================================================================\n\n";

//-----------
// Read data
//-----------

// Read Mesh data
   try {
      if (Verbosity>1)
         cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile());
      if (Verbosity > 1)
         cout << ms;
      
//    Declare problem data (boundary conditions, body forces)
      if (Verbosity>1)
         cout << "Allocating memory for matrix and R.H.S. ...\n";
      Vect<double> b(ms), u(ms);
      u.setName("Temperature");
      Prescription pr(ms,data.getDataFile());
      pr.get(INITIAL_FIELD,u);
      Vect<double> bc(ms), body_f(ms), bound_f(ms);

//    Instantiating equation
      if (Verbosity>0)
         cout << "Setting equation" << endl;
      DC2DT3 eq(ms);
      eq.setInput(INITIAL,u);

//    Loop over time steps
//    --------------------
      TimeLoop {
      
         if (Verbosity>2)
            cout << "Performing time step " << theStep << endl;

//       Read boundary temperature and sources
         pr.get(BOUNDARY_CONDITION,bc,theTime);
         pr.get(SOURCE,body_f,theTime);
         pr.get(FLUX,bound_f,theTime);

//       Prescribe boundary conditions, sources, fluxes
         eq.setInput(BOUNDARY_CONDITION,bc);
         eq.setInput(SOURCE,body_f);
         eq.setInput(FLUX,bound_f);
         eq.setTerms(LUMPED_CAPACITY|DIFFUSION);

//       Run one time step
         eq.run(TRANSIENT_ONE_STEP);

//       Save solution in plot file if necessary
         u.setTime(theTime);
         if (Verbosity>4)
            cout << "\nSolution at time: " << theTime << endl << u;
         if (theStep%save_flag == 0)
            pf.put(u);

      }

//    Calculate error for the demo
      void error(double time, const Mesh &ms, const Vect<double> &u);
      error(theTime,ms,u);
   } CATCH_EXCEPTION

   return 0;
}
