/*==============================================================================

                             *************************
                             *     C O N T A C T     *
                             *************************


                            A Signorini Contact Problem
                      for Plane Strain Linearized Elasticity

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

  ==============================================================================*/

#include "OFELI.h"
#include "Solid.h"

using namespace OFELI;

int main(int argc, char *argv[])
{
   ifstream mf, bcf, bodyf, boundf;
   ofstream dmf;
   double   penal = 1.e07;

   cout << "\n\n";
   cout << "contact, version 1.1, Copyright (c) 1998 - 2008 by Rachid Touzani\n";
   cout << "contact comes with ABSOLUTELY NO WARRANTY.\n";
   cout << "This is free software, and your are allowed to redistribute it\n";
   cout << "under certain conditions. Details are distributed with the software." << endl;

   if (argc < 2) {
      cout << "\nUsage:  contact  <parameter_file>\n";
      return 0;
   }
   IPF data("contact - 1.0",argv[1]);
   int verbose = data.getVerbose();

   if (verbose) {
      cout << endl << endl;
      cout << "    *********************************************************\n";
      cout << "    *                     C o n t a c t                     *\n";
      cout << "    * Two-Dimensional Linearized Elastostatics with contact *\n";
      cout << "    *********************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "         A Finite Element Code for Linearized Elastostatics\n";
      cout << "                      in 2-D Geometries with contact\n\n";
      cout << "          Contact uses OFELI Library of Finite Element Classes\n\n";
      cout << "                         V E R S I O N   1.1\n\n";
      cout << "                     Copyright R. Touzani, 2003\n\n";
      cout << "=====================================================================\n\n";
   }

//----------
// Read data
//----------

// Read parameters and mesh data
   try {
      int output_flag = data.getOutput();
      int save_flag = data.getSave();
      MaxNbIterations = data.getNbIter();
      theTolerance = data.getTolerance();
      Mesh ms(data.getMeshFile(1));

      Prescription p(ms,data.getDataFile());
      if (verbose > 1)
         cout << "Reading mesh data ...\n";
      if (verbose > 1)
         cout << ms;

//    Declare problem data (matrix, rhs, boundary conditions, body forces)
      if (verbose > 1)
         cout << "Allocating memory for matrix and R.H.S. ...\n";
      SkSMatrix<double> A(ms);
      Vect<double> u(ms), v(ms);

//    Read boundary conditions, body and boundary forces
      if (verbose > 1)
         cout << "Reading boundary conditions, body and boundary forces ...\n";
      Vect<double> bc(ms), body_f(ms), bound_f(ms);
      p.get(BOUNDARY_CONDITION,bc,0);
      p.get(BODY_FORCE,body_f,0);
      p.get(BOUNDARY_FORCE,bound_f,0);

//    ---------------
//    Iteration Loop
//    ---------------

      Converged = false;
      IterationLoop {
         v = 0.; A = 0.;
         MeshElements(ms) {
            Elas2DT3 eq(theElement);
            eq.Deviator();
            eq.Dilatation();
            eq.BodyRHS(body_f);
            eq.ElementAssembly(A);
            eq.ElementAssembly(v);
         }
         MeshSides(ms) {
            Elas2DT3 eq(theSide,u);
            eq.BoundaryRHS(bound_f);
            eq.SignoriniContact(bound_f,penal);
            eq.SideAssembly(A);
            eq.SideAssembly(v);
         }
         A.Prescribe(v,bc);
         A.solve(v);

         theDiscrepancy = Discrepancy(u,v,2);
         if (theDiscrepancy < theTolerance)
            Converged = true;
         if (verbose > 0)
            cout << "Iteration: " << theIteration << ", Discrepancy: " << theDiscrepancy << endl;
      }

      if (output_flag > 1)
         cout << u;
      if (save_flag) {
         IOField pl_file(data.getMeshFile(),data.getString("output_file"),ms,IOField::OUT);
         pl_file.put(u);
      }
   } CATCH_EXCEPTION
   return 0;
}
