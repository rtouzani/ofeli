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

#include "UpdateMF.h"
using namespace OFELI;

void UpdateMF(Mesh &ms0, Mesh &ms1, double omega, std::complex<double> &volt,
              Vect<std::complex<double> > &b0, Vect<std::complex<double> > &b1)
//------------------------------------------------------------------------------
//      Update magnetic field by taking into account boundary conditions
//------------------------------------------------------------------------------
{
   std::complex<double> a(0,omega);
   size_t i;

   std::complex<double> yy;
   MeshElements(ms0) {
      EC2D1T3 eq(theElement);
      yy += eq.IntegND(b0);
   }
   std::complex<double> xx;
   MeshElements(ms1) {
      EC2D1T3 eq(theElement,b1);
      xx += eq.IntegMF();
   }

// Calculate area of vacuum
   double aa=0.;
   MeshSides(ms0) {
      EC2D1T3 eq(theSide);
      aa += eq.VacuumArea();
   }
   MeshSides(ms1) {
      EC2D1T3 eq(theSide);
      aa -= eq.VacuumArea();
   }

   for (i=0; i<b0.size(); i++)
      b0[i] *= volt/(a*(xx+MU0*aa)+yy);
   for (i=0; i<b1.size(); i++)
      b1[i] *= volt/(a*(xx+MU0*aa)+yy);
}
