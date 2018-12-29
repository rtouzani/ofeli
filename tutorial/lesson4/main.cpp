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

            Solution of a 2-D Transient Heat Equation using :
            - P1 Finite elements for space discretization
            - Implicit Euler scheme for time discretization

  ==============================================================================*/


#include "OFELI.h"
#include "Therm.h"
using namespace OFELI;


int main(int argc, char *argv[])
{
// Read and output mesh data
   if (argc <= 1) {
      cout << " Usage: lesson4 <parameter_file>" << endl;
      exit(1);
   }

   try {
      IPF data(argv[1]);
      theFinalTime = data.getMaxTime();
      theTimeStep = data.getTimeStep();
      Mesh ms(data.getMeshFile());

//    Declare problem data (matrix, rhs, boundary conditions, body forces)
      SkSMatrix<double> A(ms);
      Vect<double> b(ms), u(ms), bc(ms), sf(ms);

//    Read in initial solution
      u = 0.;

//    Loop over time steps
//    --------------------
      TimeLoop {
         b = 0;
         bc = 0;
         sf.setSideBC(1,"1.0");

//       Loop over elements
         MeshElements(ms) {

//          Declare an instance of class DC2DT3
            DC2DT3 eq(theElement,u,theTime);

//          Capacity contribution to matrix and to RHS
            eq.LCapacityToLHS(1./theTimeStep);
            eq.LCapacityToRHS(1./theTimeStep);

//          Diffusion contribution to matrix
            eq.Diffusion();

//          Assemble element matrix and RHS
            if (theStep==1)
               eq.ElementAssembly(A);
            eq.ElementAssembly(b);
         }

//       Loop over edges (sides)
         MeshSides(ms) {
            DC2DT3 eq(theSide,u,theTime);
            eq.BoundaryRHS(sf);
            eq.SideAssembly(b);
         }

//       Impose boundary conditions on matrix and RHS
//       Solve the linear system of equations.
//       Factorize at first time step only
         A.Prescribe(b,bc,theStep-1);
         if (theStep==1)
            A.Factor();
         A.solve(b);
         u = b;

//       Output solution
         cout << "\nSolution for time: " << theTime << endl << u;
      }
   } CATCH_EXCEPTION

   return 0;
}
