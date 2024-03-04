/*==============================================================================

                                *******************
                                *    S T D C 2    *
                                *******************


                        A Finite Element Code for Steady-State
                Analysis of Thermal Diffusion - Convection Problems
                                  in 2-D Geometries

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
#include "Therm.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
// Expand arguments
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <param_file>\n";
      return EXIT_FAILURE;
   }
   IPF data("stdc2 - 2.0",argv[1]); 
   cout << endl << endl;
   cout << "    *******************************************************\n";
   cout << "    *                      S  T  D  C  2                  *\n";
   cout << "    *     Steady State Thermal Diffusion - Convection     *\n";
   cout << "    *                 in 2-D Geometries                   *\n";
   cout << "    *******************************************************\n\n\n";
   cout << "=====================================================================\n\n";
   cout << "               A Finite Element Code for Steady-State\n";
   cout << "           Analysis of Thermal Diffusion - Convection Problems\n";
   cout << "                             in 2-D Geometries\n\n";
   cout << "            STDC2 uses OFELI Library of Finite Element Classes\n\n";
   cout << "                           V E R S I O N   2.0\n\n";
   cout << "                     Copyright R. Touzani, 2008\n\n";
   cout << "=====================================================================\n\n";

// Read Mesh data and options
   try {
      if (Verbosity>1)
         cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile(),true);
      Prescription p(ms,data.getPrescriptionFile());
      if (Verbosity>4)
         cout << ms;

//    Declare problem data (matrix, rhs, boundary conditions, body forces)
      if (Verbosity>1)
         cout << "Allocating memory for matrix and R.H.S. ...\n";
      Vect<double> u(ms);
      u.setName("Temperature");
      if (Verbosity>1)
         cout << "Reading boundary conditions, body and boundary forces ...\n";
      Vect<double> bc(ms), body_f(ms), bound_f(ms);
      p.get(EType::BOUNDARY_CONDITION,bc);
      p.get(EType::SOURCE,body_f);
      p.get(EType::FLUX,bound_f);

//    Read velocity for convection
      Vect<double> v(ms,NODE_DOF,2);
      if (data.getInteger("v_flag")) {
         if (Verbosity>1)
            cout << "Reading Velocity in file ...\n";
         IOField ff(data.getMeshFile(),data.getString("v_file"),ms,IOField::IN);
         ff.get(v);
      }

//    Set equation features and choose solver
      DC2DT3 eq(ms,u);
      eq.setBoundaryCondition(bc);
      eq.setBodyForce(body_f);
      eq.setFlux(bound_f);
      eq.setInput(EType::VELOCITY,v);
      eq.setTerms(PDE_Terms::DIFFUSION);
      eq.setTerms(PDE_Terms::CONVECTION);
      eq.setSolver(GMRES_SOLVER,DILU_PREC);
      eq.run();

//    Output and save result
      if (data.getOutput() > 0)
         cout << u;
      if (data.getPlot())
         saveField(u,data.getPlotFile(),GMSH);
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}