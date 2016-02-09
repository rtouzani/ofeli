/*==============================================================================

                                *******************
                                *     S T D 3     *
                                *******************


                        A Finite Element Code for Steady-State
                        Analysis of Thermal Diffusion Problems
                                  in 3-D Geometries


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
#include "Therm.h"
using namespace OFELI;


int main(int argc, char *argv[])
{
// Expand arguments
   if (argc < 2) {
      cout << "\nUsage:  std3  <param_file>\n";
      return 0;
   }

   IPF data("std3 - 1.0",argv[1]);
   int output_flag = data.getOutput();
   int save_flag = data.getSave();

   if (output_flag) {
      cout << endl << endl;
      cout << "    *******************************************************\n";
      cout << "    *                         S  T  D  3                  *\n";
      cout << "    *               Steady State Thermal Diffusion     *\n";
      cout << "    *                 in 3-D Geometries                   *\n";
      cout << "    *******************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "               A Finite Element Code for Steady-State\n";
      cout << "                Analysis of Thermal Diffusion Problems\n";
      cout << "                             in 3-D Geometries\n\n";
      cout << "            STD3 uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   1.0\n\n";
      cout << "                     Copyright R. Touzani, 2002\n\n";
      cout << "=====================================================================\n\n";
   }

//---------------------------------
// Read data
//---------------------------------

// Read Mesh data
   if (output_flag > 1)
      cout << "Reading mesh data ...\n";
   Mesh ms(data.getMeshFile(),true);
   Prescription p(ms,data.getDataFile());
   if (output_flag > 1)
      cout << ms;

// Declare problem data (matrix, rhs, boundary conditions, body forces)
   if (output_flag > 1)
      cout << "Allocating memory for matrix and R.H.S. ...\n";
   SpMatrix<double> A(ms);
   Vect<double> u(ms.getNbEq()), b(ms.getNbEq());

// Read boundary conditions, body and boundary forces
   if (output_flag > 1)
      cout << "Reading boundary conditions, body and boundary sources ...\n";
   Vect<double> bc(ms), body_f(ms), bound_f(ms,1,SIDE_FIELD);
   p.get(BOUNDARY_CONDITION,bc);
   p.get(BODY_FORCE,body_f);
   p.get(BOUNDARY_FORCE,bound_f,0);

// Loop over elements
   if (output_flag > 1)
      cout << "Looping over elements ...\n";
   MeshElements(ms) {
      DC3DT4 eq(theElement);
      eq.Diffusion();
      eq.BodyRHS(body_f);
      eq.updateBC(bc);
      eq.ElementAssembly(A);
      eq.ElementAssembly(b);
   }

// Loop over sides
   if (output_flag > 1)
     cout << "Looping over sides ...\n";
   MeshSides(ms) {
      DC3DT4 eq(theSide);
      eq.BoundaryRHS(bound_f);
      eq.SideAssembly(b);
   }

   if (output_flag > 1)
      cout << "Solving linear system ...\n";
   LinearSolver<double> ls(1000,1.e-8,0);
   int nb_it = ls.solve(A,b,u,CG_SOLVER,SSOR_PREC);
   cout << "Number of iterations: " << nb_it << endl;

// Output solution
   Vect<double> v(ms);
   v.insertBC(ms,u,bc);
   if (output_flag > 0)
      cout << u;

   if (save_flag) {
      IOField pf(data.getSaveFile(),IOField::OUT);
      pf.put(v);
   }
   return 0;
}
