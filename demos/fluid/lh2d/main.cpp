/*==============================================================================

                             *************************
                             *        L H 2 D        *
                             *************************

                A Finite Volume Code for Linear Hyperbolic equations in
                                  Two Dimensions

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

  ==============================================================================
                Author: Stephane Clain
  ==============================================================================*/

#include "OFELI.h"
#include "CL.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   void set_init(Vect<double> &u);
   if (argc == 1) {
      cout << "Usage: " << argv[0] << " <parameter file>" << endl;
      return EXIT_FAILURE;
   }

   try {
      IPF proj("lh2d - 1.0",argv[1]);
      string project = proj.getProject();
      theFinalTime = proj.getMaxTime();
      int time_mod = proj.getPlot();
      Mesh ms(proj.getMeshFile());
      ms.getAllSides();
      Vect<double> u(ms,ELEMENT_DOF,1);
      LCL2DT eq(ms,u);
      Reconstruction pp(ms);

//    Get velocity field
      Vect<double> v(ms,SIDE_DOF,2);
      IOField vf(proj.getMeshFile(),proj.getString("vf"),ms,IOField::IN);
      vf.get(v);
      eq.setVelocity(v);

//    Initial condition
      set_init(u);

//    open output file
      IOField ff(proj.getMeshFile(),proj.getPlotFile(),ms,IOField::OUT);
      Vect<double> U(ms,NODE_DOF,"T",1,0.);
      pp.P0toP1(u,U);
      ff.put(U);

//    Choose MUSCL method
      eq.setMethod(Muscl::MULTI_SLOPE_M_METHOD);

      eq.setBC(proj.getDouble("bc"));
      eq.setCFL(proj.getDouble("CFL"));

//    TIME LOOP
      TimeLoop {
         eq.setReconstruction();
         eq.setBC(1,0.);                // reflection condition for the boundary side
         theTimeStep = eq.runOneTimeStep();
         if (theStep%time_mod==0) {
            cout << "Saving step " << theStep << " at time " << theTime+theTimeStep << endl;
            U.setTime(theTime+theTimeStep);
            pp.P0toP1(u,U);
            ff.put(U);
         }
      }
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}

void set_init(Vect<double> &u)
{
   ElementLoop(u.getMesh()) {
      Point<double> a = Triang3(theElement).getCenter();
      double r = sqrt((a.x-0.35)*(a.x-0.35)+(a.y-0.35)*(a.y-0.35));
      u(theElementLabel) = 0.;
      if (r<0.15)
         u(theElementLabel) = 1.;
   }
}
