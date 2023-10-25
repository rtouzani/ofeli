/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani
 
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

                 An example of a Finite Element Code using OFELI

            Solution of a 2-D parabolic problem using P1 Finite elements

  ==============================================================================*/

#include "OFELI.h"
#include "LinearPDE.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
// Expand arguments
   if (argc < 2) {
      cout << "\nUsage: " << argv[0] << " <project_file>\n";
      return EXIT_FAILURE;
   }
   cout << "A PROGRAM TO ILLUSTRATE A 2-D PARABOLIC EQUATION" << endl << endl;
   IPF data(string(argv[0])+" - 1.0",argv[1]); 

   try {

      int save_flag = data.getSave();
      Verbosity = data.getVerbose();
      theFinalTime = data.getMaxTime();
      theTimeStep = data.getTimeStep();

      Mesh ms(data.getMeshFile());
      Vect<double> u(ms);
      IOField pf("aux.pl",IOField::OUT);

//    Declare and initialize used vector
//    u: initial solution and solution at each time step
      u.set("0.");
 
//    Instantiate equation class and declare used terms
      LinearPDE2D eq(ms);

//    Build the differential system
      TimeStepping ts(BACKWARD_EULER,theTimeStep,theFinalTime);
      ts.setPDE(eq);
      ts.setInitial(u);

//    Set solver of the linear system (See class LinearSolver for other choices)
      ts.setLinearSolver(CG_SOLVER,DILU_PREC);

//    Time loop
      TimeLoop {

//       Set right-hand side of equation and boundary condition
         ts.setRHS("tanh(10*y)*(exp(t)+200*(exp(t)-1)/(cosh(10*y)^2))");
         ts.setBC(1,"tanh(10*y)*(exp(t)-1)");

//       Set pde terms for the heat equation
         eq.set_10();
         eq.set_02();

//       Run the time step: The solution is stored in vector u
         ts.runOneTimeStep();

//       Save solution
         u.setTime(theTime);
         pf.put(u);
      }

//    Save plot file to a GMSH file for plotting
//    First close file
      pf.close();
      saveGmsh(ms,"aux.pl",data.getPlotFile(),save_flag);
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}