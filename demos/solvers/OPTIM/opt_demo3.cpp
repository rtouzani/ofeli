/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2023 Rachid Touzani

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
#include "Opt3.h"
using namespace OFELI;

int main(int argc, char *argv[])
{
   if (argc <= 1) {
      cout << "Usage: " << argv[0] << " <project_file>\n";
      exit(1);
   }

   try {
      IPF proj(argv[0],argv[1]);
      Mesh ms(proj.getMeshFile(1));

//    Declare solution and boundary condition vectors
      Vect<double> u(ms), bc(ms), f(ms);
      bc.setNodeBC(1,"sin(pi*x)*cos(pi*y)");
      f = "2*pi^2*sin(pi*x)*cos(pi*y)";
 
//    Define the optimization problem and choose the optimization algorithm
      Opt3 opt(ms);
      opt.getEquation()->setInput(SOURCE,f);
      OptSolver os(opt,u);
      os.setOptMethod(OptSolver::GRADIENT);

//    Set Dirichlet boundary conditions as constraints for the optimization problem
      os.setBC(bc);

//    Run the optimization procedure
//    Some parameters are retrieved from project file
      double toler = proj.getTolerance();
      int max_it = proj.getNbIter();
      os.run(toler,max_it);

//    Output class information and error
      cout << os;
      Vect<double> sol(ms);
      sol = "sin(pi*x)*cos(pi*y)";
      cout << "L2-Norm of solution: " << (u-sol).getWNorm2() << endl;
   } CATCH_EXCEPTION

   return EXIT_SUCCESS;
}
