/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2012 Rachid Touzani

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
using namespace OFELI;

void error(double time, const Mesh &ms, const Vect<double> &u)
{
   double e=0, a=10;
   MeshNodes(ms) {
      double U = tanh(a*theNode->getCoord(2))*(exp(time)-1);
      double ee = U - u(theNodeLabel);
      e += ee*ee;
   }
   cout << "Discrete L2 error = " << sqrt(e/ms.getNbNodes()) << endl;
}
