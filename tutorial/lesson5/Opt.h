/*==============================================================================

                                    O  F  E  L  I

                            Object  Finite  Element  Library

  ==============================================================================

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

#ifndef __OPT_H
#define __OPT_H

#include "OFELI.h"
#include "Therm.h"
using namespace OFELI;

class Opt {

 public:
   Opt(Mesh &ms) { _ms = &ms; }
   void Objective(Vect<double> &x, double &f, Vect<double> &g) {
        f = 0.;
        g = 0.;
        MeshElements(*_ms) {
           DC2DT3 eq(theElement);
           Vect<double> xe(theElement,x);
           f += eq.Energy(xe);
           eq.EnerGrad(xe);
           eq.ElementAssembly(g);
        }
   }

 private:
   Mesh  *_ms;
};

#endif
