/*==============================================================================

                              **********************
                              *    T  I  F  F  2   *
                              **********************

                 A Finite Element Code for Transient Incompressible
                       Fluid Flow Simulations in 2-D Geometries

  ------------------------------------------------------------------------------

   Copyright (C) 1998 - 2016 Rachid Touzani

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
#include "User.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
      cout << "\nUsage:  tiff2  <parameter_file>\n";
      return 0;
   }

   IPF data("tiff2 - 1.0",argv[1]);
   int output_flag = data.getOutput();
   int save_flag = data.getSave();
   theTimeStep = data.getTimeStep();
   theFinalTime = data.getMaxTime();

   if (output_flag) {
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
   if (output_flag > 1)
     cout << "Reading mesh data ...\n";
   Mesh ms(data.getMeshFile());
   User ud(ms);

// Print Mesh data
   if (output_flag > 1)
      cout << ms;

// Declare problem data (matrix, rhs, boundary conditions, body forces)
   if (output_flag > 1)
      cout << "Allocating memory for matrix and R.H.S. ...\n";
   SkSMatrix<double> A(ms);

// Read initial condition, boundary conditions, body and boundary forces
   Vect<double> b(ms), u(ms), bc(ms), body_f(ms), bound_f(ms);
   Vect<double> pm(ms,1,NODE_DOF), ep(ms,1,ELEMENT_DOF), p(ms,1,NODE_DOF);
   u = 0; bc = 0; body_f = 0; bound_f = 0; 
   u.setName("Velocity");
   ud.setInitialData(u);
   ud.setDBC(bc);
   ud.setBodyForce(body_f);
   ud.setSurfaceForce(bound_f);

   IOField v_file(data.getMeshFile(),data.getString("v_file"),ms,IOField::OUT);
   IOField p_file(data.getMeshFile(),data.getString("p_file"),ms,IOField::OUT);

// Loop over time steps
// --------------------

   TimeLoop {

      if (output_flag > 1)
         cout << "Performing time step " << theStep << " ..." << endl;
      cout << "Step: " << theStep << ", Time: " << theTime << endl;
      b = 0;

//    Loop over elements
      if (output_flag > 1)
         cout << "Looping over elements ...\n";
      MeshElements(ms) {
         NSP2DQ41 eq(theElement,u,theTime);
         eq.LMass(1./theTimeStep);
         eq.Penal(1.e07);
         eq.Viscous(0.1);
         eq.RHS_Convection();
         eq.BodyRHS(ud);
         if (theStep==1)
            eq.ElementAssembly(A);
         eq.ElementAssembly(b);
      }

//    Loop over sides
      if (output_flag > 1)
        cout << "Looping over sides ...\n";
      MeshSides(ms) {
         NSP2DQ41 eq(theSide);
         eq.BoundaryRHS(ud);
         eq.SideAssembly(b);
      }

//    Impose boundary conditions and solve the linear system
      A.Prescribe(b,bc,theStep-1);
      if (theStep == 1)
         A.Factor();
      A.Solve(b);
      u = b;

//    Output and/or Store solution
      u.setTime(theTime);
      if (output_flag > 0)
         cout << u;
      if (save_flag)
         v_file.put(u);

//    Calculate and smooth pressure
      double pres = 0.;
      p.setTime(theTime);
      p.setName("Pressure");
      Reconstruction rr(ms);
      ep = 0;
      MeshElements(ms) {
         NSP2DQ41 eq(theElement,u,theTime);
         pres += eq.Pressure(1.e07);
         ep(theElementLabel) += pres;
      }
      rr.P0toP1(ep,p);
      if (output_flag > 0)
         cout << p;
      if (save_flag)
      p_file.put(p);
   }
   return 0;
}
