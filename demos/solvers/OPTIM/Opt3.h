/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

   Copyright (C) 1998 - 2021 Rachid Touzani

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

          User defined class to define the objective for an optimization
                         problem associated to a PDE

  ==============================================================================*/

#ifndef __OPT3_H
#define __OPT3_H

#include "OFELI.h"
#include "solvers/MyOpt.h"
#include "Therm.h"

using namespace OFELI;

class Opt3 : public MyOpt
{

 public:

/* Constructor using the Mesh instance
 * mesh (in) Mesh instance
 */
   Opt3(Mesh& mesh) : MyOpt(mesh)
   {
      setEquation(&eq);
      eq.setMesh(*_theMesh);
   }

/* Function to define the objective
 * x (in) Vector of optimization variables
 * Return value: Objective
 */
   double Objective(Vect<double>& x)
   {
      return eq.Energy(x);
   }

/* Function to define the gradient
 * x (in) Vector of optimization variables
 * g (out) Vector of gradient
 */
   void Gradient(Vect<double>& x,
                 Vect<double>& g)
   {
      eq.EnergyGrad(x,g);
   }

   DC2DT3 eq;
};

#endif
