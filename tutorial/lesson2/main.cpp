/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2017 Rachid Touzani

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

            Solution of a 2-D Poisson problem using P1 Finite elements

  ==============================================================================*/

#include "OFELI.h"
#include "Therm.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   banner("lesson2");
   
// Read and output mesh data
   if (argc <= 1) {
      cout << "Usage: lesson2 <mesh_file>" << endl;
      exit(1);
   }
   Mesh ms(argv[1]);

// Declare problem data (matrix, rhs, boundary conditions, body forces)
   SkSMatrix<double> A(ms);
   Vect<double> b(ms), bc(ms);
   bc.setNodeBC(1,"y");

// Loop over elements
// ------------------

   MeshElements(ms) {

//    Declare an instance of class DC2DT3 (contains equation)
      DC2DT3 eq(theElement);

//    Diffusion contribution to matrix
      eq.Diffusion();

//    Assemble element matrix and RHS
//    For this example, assembling the RHS is useless
      eq.ElementAssembly(A);
      eq.ElementAssembly(b);
   }

// Impose boundary conditions on matrix and RHS
   A.Prescribe(b,bc);

// Solve the linear system of equations
   A.solve(b);

// Output solution
   cout << "\nSolution:\n" << b;

   return 0;
}
