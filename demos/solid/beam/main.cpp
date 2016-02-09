/*==============================================================================

                                      O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

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

  ==============================================================================

                   An example of a Finite Element Code using OFELI

         Solution of a planar elastic beam problem using P1 Finite elements  

  ==============================================================================*/

#include "OFELI.h"
#include "Solid.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
     cout << "\nUsage:  beam  <parameter_file>\n";
     return 0;
   }

   IPF data("beam 1.0",argv[1]);
   int verbose = data.getVerbose();
   int output_flag = data.getOutput();

   if (output_flag) {
     cout << endl << endl;
     cout << "    *******************************************************\n";
     cout << "    *                      B  E  A  M                     *\n";
     cout << "    *                Three-Dimensional Beams              *\n";
     cout << "    *******************************************************\n\n\n";
     cout << "=====================================================================\n\n";
     cout << "            A Finite Element Code for spatial beams\n";
     cout << "           beam uses OFELI Library of Finite Element Classes\n\n";
     cout << "                        V E R S I O N   1.0\n\n";
     cout << "                    Copyright R. Touzani, 2005\n\n";
     cout << "=====================================================================\n\n";
   }

// Read Mesh data
   if (output_flag > 1)
     cout << "Reading mesh data ..." << endl;
   Mesh ms(data.getMeshFile());
   if (output_flag > 1)
      cout << ms;
   Prescription p(ms,data.getDataFile());

// Declare problem data (matrix, rhs, boundary conditions, body forces)
   SkSMatrix<double> A(ms);
   if (verbose > 1)
      cout << "Reading boundary conditions, body and boundary forces ...\n";
   Vect<double> b(ms), bc(ms), bf(ms);
   p.get(BOUNDARY_CONDITION,bc,0);
   p.get(POINT_FORCE,bf);

// Build matrix and R.H.S.
   MeshElements(ms) {
      Beam3DL2 eq(theElement,0.1,0.1,0.1);
      eq.Stiffness();
      eq.ElementAssembly(A);
   }
   b = bf;
   A.Prescribe(b,bc);
   A.solve(b);

// Output solution
   if (output_flag > 0)
      cout << b;
   Vect<double> u(ms,3);
   Beam3DL2 eq(ms,b,u);
   ms.setVerbose(10);

   cout << "Element     Axial Force       Bending Moment    Twisting Moment   Shear Force" << endl;
   MeshElements(ms) {
      Beam3DL2 eq(theElement,0.1,0.1,0.1,b);
      cout << setw(8) << theElementLabel;
      cout << " " << eq.AxialForce();
      cout << " " << eq.BendingMoment();
      cout << " " << eq.TwistingMoment();
      cout << " " << eq.ShearForce() << endl;
   }

// Transform mesh to deformed one
   DeformMesh(ms,u);
   ms.put("deformed_beam.m");
   saveTecplot("beam_tecplot.dat",ms);
   return 0;
}
