/*==============================================================================

                              **********************
                              *    T  I  N  S  2   *
                              **********************

                 A Finite Element Code for Transient Incompressible
                         Fluid Flow Simulations in 2-D

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
#include "Fluid.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
      cout << " Usage: tins2 <project_file>" << endl;
      exit(1);
   }

   try {
      IPF proj(argv[1]);
      Mesh mesh(proj.getMeshFile());
      theTimeStep = proj.getTimeStep();
      theFinalTime = proj.getMaxTime();
      int plot_flag=proj.getPlot();
      double Re=proj.getDouble("Reynolds");
      cout << endl;
      cout << "=================================================================================\n\n";
      cout << "                                   tins2\n\n";
      cout << "                       A Finite Element Code for\n";
      cout << "             2-D Dynamic Incompressible Fluid Flow Analysis\n\n";
      cout << "                          V E R S I O N   2.0\n\n";
      cout << "=================================================================================\n\n";

//    Read initial condition, boundary conditions, body and boundary forces
      IOField vff(proj.getMeshFile(),proj.getString("v_file"),mesh,IOField::OUT);
      IOField pff(proj.getMeshFile(),proj.getString("p_file"),mesh,IOField::OUT);

      Vect<double> u(mesh,"Velocity",0), p(mesh,"Pressure",0,1);
      Vect<double> bc(mesh), bf(mesh), sf(mesh);
      TINS2DT3S eq(mesh,u,p,theTimeStep,Re);
      eq.setVerbose(proj.getVerbose());
      eq.setTolerance(proj.getTolerance());
      Prescription pr(mesh,proj.getDataFile());
      pr.get(INITIAL_FIELD,u);
      eq.setInput(INITIAL_FIELD,u);

//    Loop on time steps
      TimeLoop {
         if (proj.getVerbose()>0)
            cout << "\nPerforming step: " << theStep << ", time: " << theTime << endl;
         pr.get(BOUNDARY_CONDITION,bc,theTime);
         eq.setInput(BOUNDARY_CONDITION,bc);
         pr.get(BODY_FORCE,bf,theTime);
         eq.setInput(BODY_FORCE,bf);
         eq.runOneTimeStep();
         if (plot_flag>0 && theStep%plot_flag==0) {
            p.setTime(theTime); u.setTime(theTime);
            vff.put(u); pff.put(p);
         }
      }
   } CATCH_EXCEPTION
   return 0;
}
