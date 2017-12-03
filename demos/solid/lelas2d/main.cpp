/*==============================================================================

                             *************************
                             *     L E L A S 2 D     *
                             *************************

                A Finite Element Code for Linearized Elastostatics in
                                  Two Dimensions

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
#include "Solid.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
      cout << "\nUsage:  lelas2d  <parameter_file>\n";
      return 0;
   }

   IPF data("lelas2d - 1.1 ",argv[1]);
   int output_flag = data.getOutput();
   int save_flag = data.getSave();

   if (output_flag) {
      cout << endl << endl;
      cout << "    *******************************************************\n";
      cout << "    *                   L E L A S 2 D                     *\n";
      cout << "    *     Two-Dimensional Linearized Elastostatics        *\n";
      cout << "    *******************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "         A Finite Element Code for Linearized Elastostatics\n";
      cout << "                             in 2-D Geometries\n\n";
      cout << "            LELAS2D uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   1.0\n\n";
      cout << "                       Copyright R. Touzani, 1998\n\n";
      cout << "=====================================================================\n\n";
   }

//---------------------------------
// Read data
//---------------------------------

// Read Mesh data
   if (output_flag > 1)
      cout << "Reading mesh data ...\n";
   Mesh ms(data.getMeshFile(1));
   Prescription p(ms,data.getDataFile());
   if (output_flag > 1)
      cout << ms;

// Declare problem data (matrix, rhs, boundary conditions, body forces)
   if (output_flag > 1)
      cout << "Allocating memory for matrix and R.H.S. ...\n";
   SkSMatrix<double> A(ms);
   Vect<double> b(ms);

// Read boundary conditions, body and boundary forces
   if (output_flag > 1)
     cout << "Reading boundary conditions ..." << endl;
   Vect<double> bc(ms);
   p.get(BOUNDARY_CONDITION,bc,0);
   if (output_flag > 1)
     cout << "Reading body forces ..." << endl;
   Vect<double> body_f(ms);
   p.get(BODY_FORCE,body_f);
   if (output_flag > 1)
      cout << "Reading Boundary Tractions ..." << endl;
   Vect<double> bound_f(ms,2,SIDE_DOF);
   p.get(TRACTION,bound_f,0);

// Loop over elements
// ------------------

   if (output_flag > 1)
      cout << "Looping over elements ...\n";
   MeshElements(ms) {
      Elas2DQ4 eq(theElement);
      eq.Deviator();
      eq.Dilatation();
      eq.BodyRHS(body_f);
      eq.ElementAssembly(A);
      eq.ElementAssembly(b);
   }

// Loop over sides
// ---------------

   if (output_flag > 1)
     cout << "Looping over sides ...\n";
   MeshSides(ms) {
      Elas2DQ4 eq(theSide);
      eq.BoundaryRHS(bound_f);
      eq.SideAssembly(b);
   }

// Take account for boundary conditions and solve system
// -----------------------------------------------------

   if (output_flag > 1)
     cout << "Imposing boundary conditions ...\n";
   A.Prescribe(b,bc);
   A.solve(b);

   if (output_flag > 0)
      cout << b;

   if (save_flag) {
      IOField pl_file(data.getPlotFile(),IOField::OUT);
      pl_file.put(b);
   }

// Calculate principal and Von-Mises stresses
// ------------------------------------------

   if (output_flag > 1)
      cout << "Calculating stresses ...\n";
   Vect<double> st(ms,"Principal Stress",0.0,3,ELEMENT_FIELD);
   Vect<double> vm(ms,"Von-Mises Stress",0.0,3,ELEMENT_FIELD);
   MeshElements(ms) {
      LocalVect<double,3> ste;
      Elas2DQ4 eq(theElement,b);
      eq.Stress(ste,vm(theElementLabel));
      st(theElementLabel,1) = ste(1);
      st(theElementLabel,2) = ste(2);
      st(theElementLabel,3) = ste(3);
   }
   return 0;
}
