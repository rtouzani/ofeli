/*==============================================================================
 
                                 O  F  E  L  I

                        Object  Finite  Element  Library
 
 ==============================================================================

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
 
  ==============================================================================
 
                                t s _ d e m o _ 2
 
         A program to show the usage of the class TimeStepping for a 
                   second order time dependent problem

    The program solves the elastodynamics problem in 2-D using P1 finite
         elements in space and the Newmark scheme for time integration

  ==============================================================================*/

#include "OFELI.h"
#include "Solid.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
// Read and set problem data
   theFinalTime = 1.;
   if (argc<2) {
      cout << "ts_demo2 <mesh_file> [time step]" << endl;
      exit(1);
   }
   Mesh ms(argv[1],true);
   theTimeStep = 0.1;
   if (argc>2)
      theTimeStep = atof(argv[2]);

// Declare and initialize used vector
// u: initial solution and solution at each time step
   Vect<double> u(ms), v(ms), bc(ms);
   u.set("-0.1",2);

// Instantiate equation class and declare used terms
   Elas2DT3 eq(ms);
   eq.setTerms(LUMPED_MASS|DEVIATORIC|DILATATION);

// Build the differential system
// We use the Newmark scheme
   TimeStepping ts(NEWMARK,theTimeStep,theFinalTime);
   ts.setPDE(eq);
   ts.setInitial(u,v);

// Set solver of the linear system (See class LinearSolver for other choices)
   ts.setLinearSolver(CG_SOLVER,DILU_PREC);

// Time loop
   saveField(u,ms,"u-0.vtk",VTK);
   TimeLoop {

//    Set Dirichlet boundary condition
      bc.setTime(theTime);
      ts.setBC(bc);

//    Run the time step: The solution is stored in vector u
      ts.runOneTimeStep();

//    Save solution for post processing
      saveField(u,ms,"u-"+itos(theStep)+".vtk",VTK);
      Mesh dm(ms);
      DeformMesh(dm,u,1.);
      dm.save("mm-"+itos(theStep)+".m");
   }

// Display TimeStepping class information
   cout << ts << "Solution L2-Norm: " << u.getWNorm2() << endl;
   return 0;
}
