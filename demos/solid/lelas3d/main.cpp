/*==============================================================================

                             *************************
                             *      L E L A S 3 D    *
                             *************************

                A Finite Element Code for 3-D Linearized Elastostatics

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
#include "Solid.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc < 2) {
      cout << "\nUsage:  lelas3d  <parameter_file>\n";
      return 0;
   }
   IPF data("lelas3d - 1.2",argv[1]);
   int verbose=data.getVerbose(), output_flag=data.getOutput(), save_flag=data.getSave();

   if (output_flag) {
      cout << endl << endl;
      cout << "    *******************************************************\n";
      cout << "    *                      L E L A S 3 D                  *\n";
      cout << "    *     Three-Dimensional Linearized Elastostatics      *\n";
      cout << "    *******************************************************\n\n\n";
      cout << "=====================================================================\n\n";
      cout << "         A Finite Element Code for Linearized Elastostatics\n";
      cout << "                             in 3-D Geometries\n\n";
      cout << "            lelas3d uses OFELI Library of Finite Element Classes\n\n";
      cout << "                           V E R S I O N   2.0\n\n";
      cout << "                       Copyright R. Touzani, 2011\n\n";
      cout << "=====================================================================\n\n";
   }

//----------
// Read data
//----------

// Read Mesh data
   if (output_flag > 1)
      cout << "Reading mesh data ...\n";
   Mesh ms(data.getMeshFile(),true);
   if (output_flag > 1)
      cout << ms;
   Prescription p(ms,data.getDataFile());

// Declare problem data (matrix, rhs, boundary conditions, body forces)
   if (output_flag > 1)
      cout << "Allocating memory for matrix and R.H.S. ...\n";
#if defined(USE_PETSC)
   PETScWrapper<double> w(argc-1,argv);
   PETScMatrix<double> A(ms);
   PETScVect<double> b(ms.getNbEq()), u(ms.getNbEq());
#else
   SpMatrix<double> A(ms);
   Vect<double> b(ms.getNbEq()), u(ms.getNbEq());
#endif

   if (verbose > 1)
      cout << "Reading boundary conditions, body and boundary forces ...\n";
   Vect<double> bc(ms), body_f(ms), bound_f(ms,3,SIDE_DOF);
   p.get(BOUNDARY_CONDITION,bc);
   p.get(BODY_FORCE,body_f,0.);
   p.get(BOUNDARY_FORCE,bound_f,0.);

// Loop over elements
   if (output_flag>1)
      cout << "Looping over elements ...\n";
   MeshElements(ms) {
      Elas3DT4 eq(theElement);
      eq.Deviator();
      eq.Dilatation();
      eq.BodyRHS(body_f);
      eq.updateBC(bc);
      eq.ElementAssembly(A);
      eq.ElementAssembly(b);
   }
   
// Loop over sides
   if (output_flag>1)
      cout << "Looping over sides ...\n";
   MeshSides(ms) {
      Elas3DT4 eq(theSide);
      eq.BoundaryRHS(bound_f);
      eq.SideAssembly(b);
   }

// Solve system
#if defined(USE_PETSC)
   w.setLinearSystem(A,b,KSPCG,PCILU,1.e-8);
   w.solve(u);
   int nb_it = w.getIterationNumber();
#else
   LinearSolver<double> ls(1000,1.e-8,0);
   int nb_it = ls.solve(A,b,u,CG_SOLVER,SSOR_PREC);
#endif
   cout << "Number of iterations: " << nb_it << endl;

#if defined(USE_PETSC)
   PETScVect<double> uf(ms);
#else
   Vect<double> uf(ms);
#endif
   uf.insertBC(u,bc);
   if (output_flag > 0)
      cout << uf;

   if (save_flag) {
      IOField pf(data.getPlotFile(),IOField::OUT);
      pf.put(uf);
      DeformMesh(ms,uf,1.e-3);
      ms.put(data.getProject()+"-1.m");
   }
   return 0;
}
