/*==============================================================================

                             *******************************
                             *        E u l e r 2 D        *
                             *******************************

                A Finite Volume Code for Compressible Euler equations in
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

  ==============================================================================
                          Author: Stephane Clain
             Sod Tube 2D with multislope method (clain 05/2007)
 ====================================================================================*/

#include "OFELI.h"
#include "CL.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
      cout << " Usage: euler2d <project file>" << endl;
      exit(1);
   }

   try {
      IPF proj("euler - 2.0",argv[1]);

      string project = proj.getProject();
      theFinalTime = proj.getMaxTime();
      Mesh ms(proj.getMeshFile());
      ms.getAllSides();

//    Read boundary conditions
      LocalVect<double,4> L, R;
      proj.get("LeftDensity",L(1));
      proj.get("LeftXVelocity",L(2));
      proj.get("LeftYVelocity",L(3));
      proj.get("LeftPressure",L(4));
      proj.get("RightDensity",R(1));
      proj.get("RightXVelocity",R(2));
      proj.get("RightYVelocity",R(3));
      proj.get("RightPressure",R(4));
      int time_mod = proj.getPlot();

      Vect<double> r(ms,1,ELEMENT_DOF), p(ms,1,ELEMENT_DOF), v(ms,2,ELEMENT_DOF);
      ICPG2DT eq(ms,r,v,p);
      eq.setInitialConditionShockTube(L,R,0.5);

      IOField ff_r(proj.getString("r-file"),IOField::OUT),
              ff_p(proj.getString("p-file"),IOField::OUT),
              ff_v(proj.getString("v-file"),IOField::OUT),
              ff_c(proj.getString("c-file"),IOField::OUT);

      Vect<double> aux1(ms), aux2(ms,2), sound(ms,1,ELEMENT_DOF);
      Reconstruction rr(ms);

//    save variables
      theTime = 0.;
      aux1.setTime(theTime);
      aux1.setName("scalar");
      aux2.setTime(theTime);
      aux2.setName("vector");
      rr.P0toP1(r,aux1);
      ff_r.put(aux1);
      rr.P0toP1(p,aux1);
      ff_p.put(aux1);
      rr.P0toP1(v,aux2);
      ff_v.put(aux2);
      eq.getSoundSpeed(sound);
      rr.P0toP1(sound,aux1);
      ff_c.put(aux1);

//    Choose solver 
      eq.setSolver(ICPG2DT::RUSANOV_SOLVER);

//    Choose MUSCL method
      eq.setMethod(Muscl::MULTI_SLOPE_M_METHOD);
      eq.setCFL(proj.getDouble("CFL"));

      MeshBoundarySides(ms)
         theSide->setCode(1,1);

      TimeLoop {
         eq.setReconstruction();
         eq.setBC(1,1.); // reflection condition for the boundary side
         cout << "Performing step " << theStep << " at time: " 
              << theTime+theTimeStep << " ..." << endl;
         theTimeStep = eq.runOneTimeStep();
         if (!(theStep%time_mod)) {
            cout << "Saving step " << theStep << " at time: " << theTime+theTimeStep << endl;
            aux1.setTime(theTime+theTimeStep);
            aux2.setTime(theTime+theTimeStep);
            rr.P0toP1(r,aux1);
            ff_r.put(aux1);
            rr.P0toP1(p,aux1);
            ff_p.put(aux1);
            rr.P0toP1(v,aux2);
            ff_v.put(aux2);
            eq.getSoundSpeed(sound);
            rr.P0toP1(sound,aux1);
            ff_c.put(aux1);
         }
      }
   } CATCH_EXCEPTION
   return 0;
}
