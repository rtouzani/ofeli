/*==============================================================================

                             *************************
                             *      L E L A S 3 D    *
                             *************************

                A Finite Element Code for 3-D Linearized Elastostatics

  ------------------------------------------------------------------------------

   Copyright (C) 1998 - 2026 Rachid Touzani

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
      cout << "\nUsage: " << argv[0] << " <parameter_file>\n";
      return EXIT_FAILURE;
   }
   IPF data("lelas3d - 1.2",argv[1]);
   int save_flag=data.getSave();
   Verbosity = data.getVerbose();

   if (Verbosity) {
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
   try {
      if (Verbosity > 1)
         cout << "Reading mesh data ...\n";
      Mesh ms(data.getMeshFile(),true);
      if (Verbosity > 5)
         cout << ms;
      Prescription p(ms,data.getDataFile());
#if defined(USE_PETSC)
      PETScWrapper<double> w(argc-1,argv);
      PETScMatrix<double> A(ms);
      PETScVect<double> b(ms.getNbEq()), u(ms.getNbEq());
#endif

      if (Verbosity > 1)
         cout << "Reading boundary conditions, body and boundary forces ...\n";
      Vect<double> u(ms), sf(ms,3,BOUNDARY_SIDE_DOF);
      p.getTraction(sf);

//    Solve equation
      if (Verbosity > 1)
         cout << "Setting and solving equation ...\n";
      Elas3DT4 eq(ms,u);
      eq.setBoundaryCondition(p.getBoundaryCondition());
      eq.setBodyForce(p.getBodyForce());
      eq.setTraction(sf);
      eq.setSolver(CG_SOLVER,DILU_PREC);
      eq.run();

#if defined(USE_PETSC)
      w.setLinearSystem(A,b,KSPCG,PCILU,1.e-8);
      w.solve(u);
      int nb_it = w.getIterationNumber();
      cout << "Number of iterations: " << nb_it << endl;
#endif

#if defined(USE_PETSC)
      PETScVect<double> uf(ms);
#endif
      if (Verbosity > 3)
         cout << u;

      if (save_flag) {
         IOField pf(data.getPlotFile(),IOField::OUT);
         pf.put(u);
         ms.Deform(u);
         ms.put(data.getProject()+"-1.m");
      }
   } CATCH_EXCEPTION
   return EXIT_SUCCESS;
}
