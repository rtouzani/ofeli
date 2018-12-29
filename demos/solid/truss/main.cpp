/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

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

               An example of a Finite Element Code using OFELI

                  Solution of a structural truss problem

  ==============================================================================*/

#include "OFELI.h"
#include "Solid.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc) { }
   IPF data("truss",argv[1]);
   int verbose = data.getVerbose();
   int output_flag = data.getOutput();

   if (output_flag) {
      cout << endl << endl;
      cout << "    *******************************************************\n";
      cout << "    *                      T R U S S                      *\n";
      cout << "    *   Two-Dimensional Trusses in Structural Mechanics   *\n";
      cout << "    *******************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "         A Finite Element Code for planar trusses of bars\n";
      cout << "              beam uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   1.0\n\n";
      cout << "                       Copyright R. Touzani, 2008\n\n";
      cout << "=====================================================================\n\n";
   }

// Read Mesh data
   try {
      if (output_flag > 1)
         cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile());
      if (output_flag > 1)
         cout << ms;
      Prescription p(ms,data.getDataFile());

//    Declare problem data (matrix, rhs, boundary conditions, body forces)
      SkSMatrix<double> A(ms);

      if (verbose > 1)
         cout << "Reading boundary conditions, body and boundary forces ...\n";
      Vect<double> bc(ms), bf(ms);
      p.get(BOUNDARY_CONDITION,bc);
      p.get(POINT_FORCE,bf);
      Vect<double> b(bf);
      double section = data.getDouble("section");

//    Loop over elements
      if (verbose > 1)
         cout << "Calculating and assembling element arrays ..." << endl;
      MeshElements(ms) {

//       Declare an instance of class Bar2DL2
         Bar2DL2 eq(theElement,section);

//       Stiffness matrix
         eq.Stiffness();

//       Assemble element matrix
         eq.ElementAssembly(A);
      }

//    Impose a concentrated vertical load at node 2
      A.Prescribe(b);

//    Solve the linear system of equations
      if (verbose > 1)
         cout << "Solving the linear system ...\n";
      A.solve(b);

//    Output solution
      if (verbose > 1)
         cout << "\nSolution:\n" << b;

//    Transform mesh to deformed one
      IOField pf(data.getMeshFile(),data.getPlotFile(),ms,IOField::OUT);
      pf.put(b);
      DeformMesh(ms,b,40);
      ms.put(data.getProject()+"-1.m");
   } CATCH_EXCEPTION
   return 0;
}
