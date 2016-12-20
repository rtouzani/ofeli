/*==============================================================================

               A Finite Element Code to Solve Helmholtz Equation
                      in a Bounded Domain using OFELI

  ------------------------------------------------------------------------------

   Copyright (C) 1998 - 2008 Rachid Touzani

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
#include "User.h"
using namespace OFELI;

void error(const Mesh &ms, User &ud, const Vect<complex<double> > &u)
{
   double ee=0, e;
   MeshNodes(ms) {
      complex<double> U = ud.ExactSolution(theNode->getCoord());
      e = Abs(U - u(theNode->n()));
      ee += e*e;
   }
   cout << "Discrete L2 error = " << sqrt(ee/ms.getNbNodes()) << endl;
}
