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

        User defined class to define the objective for an optimization
                      problem associated to a PDE

  ==============================================================================*/

#ifndef __OPT_H
#define __OPT_H

#include "OFELI.h"
#include "solvers/MyOpt.h"
#include "Therm.h"
using namespace OFELI;

class Opt : public MyOpt
{

 public:

/* Constructor using the Mesh instance
 * mesh (in) Mesh instance
 */
   Opt(const Mesh& mesh) : MyOpt(mesh)
   { }


/* Function to define the objective
 * x (in) Vector of optimization variables
 * Return value: Objective
 */
   double Objective(const Vect<double>& x)
   {
      double f = 0.;
      MeshElements(*_theMesh) {
         DC2DT3 eq(theElement);
         f += eq.Energy(Vect<double>(theElement,x));
      }
      return f;
   }

/* Function to define the gradient
 * x (in) Vector of optimization variables
 * g (out) Vector of gradient
 */
   void Gradient(const Vect<double>& x,
                 Vect<double>&       g)
   {
      g = 0.;
//    The gradient is assembled element by element
      MeshElements(*_theMesh) {
         DC2DT3 eq(theElement);
         eq.EnerGrad(Vect<double>(theElement,x));
         eq.ElementAssembly(g);
      }
   }

};

#endif
