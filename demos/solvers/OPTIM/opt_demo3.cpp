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

               A program to illustrate the Optimization solver

                 Solution of the Laplace equation formulated
                           as an optimization problem
 
  ==============================================================================*/

#include "OFELI.h"
#include "Opt.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc <= 1) {
      cout << " Usage: opt_demo3 <project_file>\n";
      exit(1);
   }

   try {
      IPF proj(argv[1]);
      Mesh ms(proj.getMeshFile(1));
      int n = ms.getNbDOF();

//    Declare solution and boundary condition vectors
      Vect<double> x(n), bc(ms);
      bc.setNodeBC(2,"1");
 
//    Define the optimization problem and choose the optimization algorithm
      Opt opt(ms);
      OptSolver os(opt,x);
      os.setOptMethod(OptSolver::TRUNCATED_NEWTON);

//    Set Dirichlet boundary conditions as constraints for the optimization problem
      os.setBC(bc);

//    Run the optimization procedure
//    Some parameters are retrieved from project file
      double toler = proj.getTolerance();
      int max_it = proj.getNbIter();
      int verb = proj.getVerbose();
      os.run(toler,max_it,verb);

//    Output class information and solution
      cout << os;
      cout << "\nSolution:\n" << x;
   } CATCH_EXCEPTION

   return 0;
}
